// Shared includes for PALMER pipeline components.
#ifndef PALMER_SCP_COMMON_HPP
#define PALMER_SCP_COMMON_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <signal.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

using namespace std;

inline int run_process(const vector<string> &args, const string &output_path = "",
                       const vector<pair<string, string>> &env_vars = {}) {
    if (args.empty()) {
        return -1;
    }

    pid_t pid = fork();
    if (pid == 0) {
        for (const auto &env : env_vars) {
            setenv(env.first.c_str(), env.second.c_str(), 1);
        }

        if (!output_path.empty()) {
            int fd = open(output_path.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644);
            if (fd == -1) {
                _exit(127);
            }
            dup2(fd, STDOUT_FILENO);
            close(fd);
        }

        vector<char *> argv;
        argv.reserve(args.size() + 1);
        for (const auto &arg : args) {
            argv.push_back(const_cast<char *>(arg.c_str()));
        }
        argv.push_back(nullptr);

        execvp(args[0].c_str(), argv.data());
        _exit(127);
    }

    int status = 0;
    waitpid(pid, &status, 0);
    if (!output_path.empty() && status != 0) {
        remove(output_path.c_str());
    }
    return status;
}

inline bool stream_process_output(const vector<string> &args,
                                  const function<bool(const char *, ssize_t)> &consumer,
                                  const vector<pair<string, string>> &env_vars = {}) {
    if (args.empty()) {
        return false;
    }

    int pipefd[2];
    if (pipe(pipefd) != 0) {
        return false;
    }

    pid_t pid = fork();
    if (pid == 0) {
        for (const auto &env : env_vars) {
            setenv(env.first.c_str(), env.second.c_str(), 1);
        }
        dup2(pipefd[1], STDOUT_FILENO);
        close(pipefd[0]);
        close(pipefd[1]);

        vector<char *> argv;
        argv.reserve(args.size() + 1);
        for (const auto &arg : args) {
            argv.push_back(const_cast<char *>(arg.c_str()));
        }
        argv.push_back(nullptr);

        execvp(args[0].c_str(), argv.data());
        _exit(127);
    }

    close(pipefd[1]);
    char buffer[8192];
    ssize_t count = 0;
    bool keep_reading = true;
    while (keep_reading && (count = read(pipefd[0], buffer, sizeof(buffer))) > 0) {
        keep_reading = consumer(buffer, count);
    }
    close(pipefd[0]);

    int status = 0;
    if (!keep_reading) {
        kill(pid, SIGTERM);
    }
    waitpid(pid, &status, 0);
    return status == 0;
}

#endif // PALMER_SCP_COMMON_HPP
