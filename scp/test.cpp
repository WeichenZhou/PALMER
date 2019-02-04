//copyright by ArthurZhou @ UMich&FUDAN&HUST
#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cstdlib>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

using namespace std;

int main(){
    
    ifstream file1;
    file1.open("test.txt");
    
    string input;
    getline(file1,input);
    cout<<input<<endl;
    file1.close();
    file1.clear();
    file1.open("test.txt");
    getline(file1,input,'\t');
    cout<<input<<endl;
    file1.close();
    file1.clear();
    file1.open("test.txt");
    getline(file1,input,' ');
    cout<<input<<endl;
}
