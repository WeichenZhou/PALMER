"""Exercise make's HTSlib requirement without installing build dependencies.

Run with ``python3 -m unittest discover -s tests -p test_build_requirements.py -v``.
Each test copies the real makefile into a temporary directory and substitutes
pkg-config and the compiler. No PALMER binary or repository file is generated.
"""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
MAKE = shutil.which("make")


@unittest.skipUnless(MAKE, "These tests require make.")
class BuildRequirementTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="palmer-build-requirements-")
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        shutil.copyfile(REPO_ROOT / "makefile", self.directory / "makefile")
        (self.directory / "PALMER.cpp").write_text("int main() { return 0; }\n")
        self.compiler_log = self.directory / "compiler-invocations.jsonl"
        self.pkg_config = self.directory / "fake-pkg-config"
        self.compiler = self.directory / "fake-cxx"
        self.write_executable(
            self.pkg_config,
            """import os
import re
import sys

args = sys.argv[1:]
version = os.environ.get("FAKE_HTSLIB_VERSION", "")
if not version:
    sys.stderr.write("Package htslib was not found.\\n")
    sys.exit(1)

def version_parts(value):
    return tuple(int(part) for part in re.findall(r"\\d+", value))

minimum = next((arg.split("=", 1)[1] for arg in args
                if arg.startswith("--atleast-version=")), None)
if "--atleast-version" in args:
    minimum = args[args.index("--atleast-version") + 1]
if minimum is not None:
    sys.exit(0 if version_parts(version) >= version_parts(minimum) else 1)
if "--modversion" in args:
    print(version)
elif "--cflags" in args or "--libs" in args:
    flags = []
    if "--cflags" in args:
        flags.append(os.environ.get("FAKE_HTSLIB_CFLAGS", "-I/fake/htslib/include"))
    if "--libs" in args:
        flags.append("-L/fake/htslib/lib -lhts")
    print(" ".join(flags))
elif "--exists" not in args:
    sys.stderr.write("Unsupported fake pkg-config query: " + repr(args) + "\\n")
    sys.exit(2)
""",
        )
        self.write_executable(
            self.compiler,
            """import json
import os
from pathlib import Path
import sys

args = sys.argv[1:]
with open(os.environ["FAKE_COMPILER_LOG"], "a") as log:
    log.write(json.dumps(args) + "\\n")
if "-o" not in args:
    sys.exit(2)
Path(args[args.index("-o") + 1]).write_text("fake PALMER executable\\n")
""",
        )

    @staticmethod
    def write_executable(path, source):
        path.write_text("#!" + sys.executable + "\n" + source)
        path.chmod(0o755)

    def run_make(self, *targets, version="1.24", cflags=None, missing_tool=False):
        environment = os.environ.copy()
        # Parent make invocations and local dependency settings must not leak
        # into the isolated build's observable acceptance/rejection behavior.
        for name in ("MAKEFLAGS", "MFLAGS", "MAKELEVEL"):
            environment.pop(name, None)
        environment["FAKE_HTSLIB_VERSION"] = version
        environment["FAKE_COMPILER_LOG"] = str(self.compiler_log)
        environment.pop("FAKE_HTSLIB_CFLAGS", None)
        if cflags is not None:
            environment["FAKE_HTSLIB_CFLAGS"] = cflags
        pkg_config = self.directory / "missing-pkg-config" if missing_tool else self.pkg_config
        return subprocess.run(
            [MAKE, "--no-print-directory", "PKG_CONFIG=" + str(pkg_config),
             "CXX=" + str(self.compiler), *targets],
            cwd=self.directory,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
        )

    def compiler_invocations(self):
        if not self.compiler_log.exists():
            return []
        return [json.loads(line) for line in self.compiler_log.read_text().splitlines()]

    def existing_binary(self):
        target = self.directory / "PALMER"
        target.write_text("existing PALMER executable\n")
        previous = time.time() - 60
        os.utime(self.directory / "PALMER.cpp", (previous, previous))
        return target

    def assert_rejected_without_compiling(self, result):
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertEqual(self.compiler_invocations(), [], result.stdout)

    def test_requires_pkg_config_executable(self):
        result = self.run_make(missing_tool=True)
        self.assert_rejected_without_compiling(result)
        self.assertIn("pkg-config", result.stdout.lower())

    def test_requires_htslib_package(self):
        result = self.run_make(version="")
        self.assert_rejected_without_compiling(result)
        self.assertIn("htslib", result.stdout.lower())

    def test_rejects_htslib_older_than_1_24(self):
        for version in ("1.23", "1.23.2"):
            with self.subTest(version=version):
                result = self.run_make(version=version)
                self.assert_rejected_without_compiling(result)
                self.assertIn("1.24", result.stdout)
                self.assertIn(version, result.stdout)

    def test_accepts_htslib_1_24(self):
        result = self.run_make(version="1.24")
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual(len(self.compiler_invocations()), 1, result.stdout)
        self.assertTrue((self.directory / "PALMER").exists())

    def test_accepts_newer_htslib_with_empty_include_flags(self):
        result = self.run_make(version="1.25", cflags="")
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual(len(self.compiler_invocations()), 1, result.stdout)

    def test_accepts_minimum_htslib_with_empty_include_flags(self):
        result = self.run_make(version="1.24", cflags="")
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual(len(self.compiler_invocations()), 1, result.stdout)

    def test_clean_works_without_dependency_tools(self):
        for missing_tool in (False, True):
            with self.subTest(missing_tool=missing_tool):
                target = self.existing_binary()
                result = self.run_make("clean", version="", missing_tool=missing_tool)
                self.assertEqual(result.returncode, 0, result.stdout)
                self.assertFalse(target.exists())
                self.assertEqual(self.compiler_invocations(), [], result.stdout)

    def test_existing_binary_cannot_bypass_minimum_htslib_version(self):
        target = self.existing_binary()
        original = target.read_bytes()
        result = self.run_make(version="1.23")
        self.assert_rejected_without_compiling(result)
        self.assertIn("1.24", result.stdout)
        self.assertEqual(target.read_bytes(), original)

    def test_supported_dependency_does_not_force_rebuilding_existing_binary(self):
        first = self.run_make(version="1.24")
        self.assertEqual(first.returncode, 0, first.stdout)
        second = self.run_make(version="1.24")
        self.assertEqual(second.returncode, 0, second.stdout)
        self.assertEqual(len(self.compiler_invocations()), 1, second.stdout)

    def test_changed_makefile_rebuilds_existing_binary(self):
        first = self.run_make(version="1.24")
        self.assertEqual(first.returncode, 0, first.stdout)
        previous = time.time() - 60
        os.utime(self.directory / "PALMER.cpp", (previous, previous))
        os.utime(self.directory / "PALMER", (previous + 10, previous + 10))
        with (self.directory / "makefile").open("a") as makefile:
            makefile.write("\n# Updated build requirements.\n")
        second = self.run_make(version="1.24")
        self.assertEqual(second.returncode, 0, second.stdout)
        self.assertEqual(len(self.compiler_invocations()), 2, second.stdout)


if __name__ == "__main__":
    unittest.main()
