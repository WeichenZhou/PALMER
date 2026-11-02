"""Regression tests for no-TSD records and unique alignment support counts.

Run with ``python3 -m unittest discover -s tests -v``. A C++17 compiler is
required: use CXX or an available clang++, g++, or c++. HTSlib, BLAST,
reference FASTA files, and third-party Python packages are not required.
The temporary harness compiles the real stage 6
and stage 8 sources. A schema-checked Python substitute performs only stage
7's no-TSD table conversion; these tests do not exercise the actual fp_ex
binary or its reference/junction checks. Synthetic negative junction results
also exercise the TSD-on caller's no-confident-TSD fallback.

PALMER_TEST_TSD_SOURCE can select a historical stage 6 source file for a
regression check. PALMER_TEST_CALLING_SOURCE similarly selects stage 8.
Each selected file must have its common.hpp available beside it.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Candidate:
    name: str
    locus: int
    strand: str = "+"


@dataclass(frozen=True)
class SupportRow:
    """One candidate HSP; repeated alignment indices retain the same full tag."""

    alignment: int
    query_start: int = 1
    query_end: int = 300
    read_start: int = 100
    read_end: int = 400


class CustomNoTsdRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler_command = shlex.split(os.environ.get("CXX", ""))
        if not compiler_command:
            for name in ("clang++", "g++", "c++"):
                compiler = shutil.which(name)
                if compiler:
                    compiler_command = [compiler]
                    break
            if not compiler_command:
                raise RuntimeError("These tests require a C++17 compiler; set CXX.")

        cls.build_directory = tempfile.TemporaryDirectory(prefix="palmer-no-tsd-build-")
        cls.addClassCleanup(cls.build_directory.cleanup)
        build_path = Path(cls.build_directory.name)
        tsd_source = Path(
            os.environ.get(
                "PALMER_TEST_TSD_SOURCE", str(REPO_ROOT / "scp" / "6_TSD_seq.cpp")
            )
        ).expanduser().resolve()
        calling_source = Path(
            os.environ.get(
                "PALMER_TEST_CALLING_SOURCE", str(REPO_ROOT / "scp" / "8_calling.cpp")
            )
        ).expanduser().resolve()
        harness = build_path / "no_tsd_harness.cpp"
        harness.write_text(
            f"#include {json.dumps(str(tsd_source))}\n"
            f"#include {json.dumps(str(calling_source))}\n"
            "int main(int argc, char **argv) {\n"
            "    if (argc != 3) return 2;\n"
            "    const string stage = argv[1];\n"
            "    string directory = argv[2];\n"
            "    if (directory.empty() || directory.back() != '/') directory += '/';\n"
            '    if (stage == "tsd") return tsd_module(directory, "custom-mtdna.fa", 0);\n'
            '    if (stage == "tsd_on") return tsd_module(directory, "custom-mtdna.fa", 1);\n'
            '    if (stage == "calling") return calling(directory, "custom-mtdna.fa", 0);\n'
            '    if (stage == "calling_on") return calling(directory, "custom-mtdna.fa", 1);\n'
            '    if (stage == "calling_sva") return calling(directory, "SVA", 0);\n'
            '    if (stage == "calling_sva_on") return calling(directory, "SVA", 1);\n'
            "    return 2;\n"
            "}\n"
        )
        cls.executable = build_path / "no_tsd_harness"
        result = subprocess.run(
            compiler_command + [
                "-std=c++17",
                "-O0",
                str(harness),
                "-o",
                str(cls.executable),
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode:
            raise RuntimeError(f"Could not compile regression harness:\n{result.stderr}")

    def run_stage(self, stage, directory):
        result = subprocess.run(
            [str(self.executable), stage, str(directory)],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    @staticmethod
    def write_inputs(directory, candidates):
        """Make small aligned-read fixtures, with identical query spans per locus."""
        cigar_rows = []
        sam_rows = []
        candidate_rows = []
        expected_ids = []
        for ordinal, candidate in enumerate(candidates, 1):
            sam_position = 1000 + ordinal * 500
            flag = 16 if candidate.strand == "-" else 0
            tagged_name = f"{candidate.name}_{sam_position}_{ordinal}"
            cigar_rows.append(
                f"500M0E\t{ordinal}\t{candidate.name}\tchr1\t{sam_position}\t500\n"
            )
            sam_rows.append(
                f"{candidate.name}\t{flag}\tchr1\t{sam_position}\t60\t500M"
                f"\t*\t0\t0\t{'ACGT' * 125}\t*\n"
            )
            fields = [
                tagged_name,
                "1", "100",       # custom-sequence coordinates
                "100", "200",     # read coordinates
                str(candidate.locus), str(candidate.locus),
                "chr1", candidate.strand,
                "500",              # read length
                *(["0"] * 7),        # no inverted segment
            ]
            candidate_rows.append("\t".join(fields) + "\n")
            expected_ids.append(".".join(fields[:9] + fields[10:]))

        (directory / "cigar.2").write_text("".join(cigar_rows))
        (directory / "region.sam").write_text("".join(sam_rows))
        (directory / "read_result.txt").write_text("".join(candidate_rows))
        return expected_ids

    def convert_no_tsd_rows(self, directory, expected_ids):
        """Substitute only stage 7's six-coordinate no-TSD format conversion."""
        rows = [
            line.split("\t")
            for line in (directory / "TSD_blastn_pre.txt").read_text().splitlines()
        ]
        for row in rows:
            self.assertEqual(
                len(row), 7,
                "A dummy TSD record must contain its ID and exactly six numeric fields.",
            )
            self.assertEqual(row[1:], ["0"] * 6)
        self.assertEqual(
            [row[0] for row in rows], expected_ids,
            "Every candidate must retain its complete ID in the dummy TSD table.",
        )
        # fp_ex(..., tsd_index=0) copies the six coordinates and appends these
        # four fields. The actual fp_ex source requires HTSlib and is not built.
        converted = ["\t".join(row + ["0", "1", "1", "N"]) + "\n" for row in rows]
        (directory / "TSD_blastn.txt").write_text("".join(converted))

    def check_candidates(self, candidates, directory=None):
        if directory is None:
            with tempfile.TemporaryDirectory(prefix="palmer-no-tsd-case-") as temporary:
                self.check_candidates(candidates, Path(temporary))
            return

        expected_ids = self.write_inputs(directory, candidates)
        self.run_stage("tsd", directory)
        self.convert_no_tsd_rows(directory, expected_ids)
        self.run_stage("calling", directory)

        expected_counts = Counter(("chr1", c.locus, c.strand) for c in candidates)
        expected_names = defaultdict(set)
        for candidate in candidates:
            expected_names[("chr1", candidate.locus, candidate.strand)].add(candidate.name)

        calls = [line.split("\t") for line in (directory / "calls.txt").read_text().splitlines()]
        self.assertEqual(len(calls), len(expected_counts))
        cluster_keys = {}
        for row in calls:
            self.assertEqual(len(row), 24)
            key = (row[1], int(row[2]), row[16])
            self.assertIn(key, expected_counts)
            cluster_keys[row[0]] = key
            self.assertEqual(row[10], "0", "Disabling TSD finding must not create confident TSD support.")
            self.assertEqual(int(row[11]), expected_counts[key])
            self.assertEqual(int(row[12]), expected_counts[key])
            self.assertEqual(row[13:16], ["0", "0", "0"])
            self.assertEqual(row[17:21], ["-", "0", "0", "0"])
        self.assertEqual(set(cluster_keys.values()), set(expected_counts))

        listed_names = defaultdict(set)
        go_through_names = defaultdict(set)
        for line in (directory / "all_reads_output.txt").read_text().splitlines():
            row = line.split("\t")
            self.assertEqual(len(row), 5)
            self.assertIn(row[0], cluster_keys)
            key = cluster_keys[row[0]]
            listed_names[key].add(row[2])
            if row[1] == "potential_go_through":
                go_through_names[key].add(row[2])
        # The existing caller may emit a seed read twice under different
        # labels. Compare identities rather than asserting row uniqueness.
        self.assertEqual(dict(listed_names), dict(expected_names))
        self.assertEqual(dict(go_through_names), dict(expected_names))

    def test_no_candidates(self):
        self.check_candidates([])

    def test_single_candidate(self):
        self.check_candidates([Candidate("single_read", 200)])

    def test_all_candidates_in_one_cluster(self):
        self.check_candidates([Candidate(f"read_{i}", 200) for i in range(3)])

    def test_candidates_in_multiple_clusters(self):
        self.check_candidates(
            [Candidate(f"read_{i}", locus) for i, locus in enumerate([200, 200, 10000, 10000, 20000, 20000])]
        )

    def test_mixed_strands_and_input_order(self):
        candidates = [
            Candidate("forward_a", 200, "+"),
            Candidate("forward_b", 200, "+"),
            Candidate("reverse_a", 200, "-"),
            Candidate("reverse_b", 200, "-"),
            Candidate("distant_a", 10000, "+"),
            Candidate("distant_b", 10000, "+"),
        ]
        orders = [
            candidates,
            list(reversed(candidates)),
            [candidates[i] for i in [2, 4, 0, 5, 1, 3]],
        ]
        for ordinal, ordered in enumerate(orders):
            with self.subTest(order=ordinal):
                self.check_candidates(ordered)

    def test_reused_regional_directory_discards_previous_records(self):
        with tempfile.TemporaryDirectory(prefix="palmer-no-tsd-reused-") as temporary:
            directory = Path(temporary)
            runs = [
                [Candidate(f"original_{i}", 200) for i in range(3)],
                [Candidate("replacement_read", 10000, "-")],
                [],
            ]
            for ordinal, candidates in enumerate(runs):
                with self.subTest(run=ordinal):
                    self.check_candidates(candidates, directory)

    def test_tsd_enabled_with_no_candidates(self):
        with tempfile.TemporaryDirectory(prefix="palmer-tsd-on-empty-") as temporary:
            directory = Path(temporary)
            self.write_inputs(directory, [])
            (directory / "TSD_blastn_pre.txt").write_text("stale_previous_record\n")
            self.run_stage("tsd_on", directory)
            for filename in (
                "TSD_blastn_pre.txt", "read_result_TSD.txt",
                "read_result_ins_seq.txt", "read_result_junction.txt",
            ):
                self.assertEqual((directory / filename).read_text(), "")

    def support_counts(self, support_rows, tsd_enabled=False, sva=False, keep_tsd_ids=True):
        """Run real stages 6/8 with repeated HSPs from explicitly shared alignments."""
        with tempfile.TemporaryDirectory(prefix="palmer-support-counts-") as temporary:
            directory = Path(temporary)
            alignments = [
                Candidate(f"alignment_{i}", 200)
                for i in range(max(row.alignment for row in support_rows) + 1)
            ]
            self.write_inputs(directory, alignments)
            original_rows = [
                line.split("\t")
                for line in (directory / "read_result.txt").read_text().splitlines()
            ]
            rows = []
            for support in support_rows:
                row = original_rows[support.alignment].copy()
                row[1:5] = map(str, (
                    support.query_start, support.query_end,
                    support.read_start, support.read_end,
                ))
                rows.append(row)
            (directory / "read_result.txt").write_text(
                "".join("\t".join(row) + "\n" for row in rows)
            )
            expected_ids = [".".join(row[:9] + row[10:]) for row in rows]
            self.run_stage("tsd", directory)
            self.convert_no_tsd_rows(directory, expected_ids)

            if tsd_enabled:
                # Synthetic caller input: negative junction checks force the
                # TSD-on caller's no-confident-TSD fallback. This does not run
                # or simulate actual BLAST-based TSD/FP validation.
                fp_rows = [
                    line.split("\t")
                    for line in (directory / "TSD_blastn.txt").read_text().splitlines()
                ]
                for row in fp_rows:
                    row[8:10] = ["0", "0"]
                (directory / "TSD_blastn.txt").write_text(
                    "".join("\t".join(row) + "\n" for row in fp_rows)
                )
            if not keep_tsd_ids:
                (directory / "TSD_blastn.txt").write_text("")
            stage = "calling_sva" if sva else "calling"
            self.run_stage(stage + ("_on" if tsd_enabled else ""), directory)
            calls = [
                line.split("\t") for line in (directory / "calls.txt").read_text().splitlines()
            ]
            self.assertTrue(calls)
            for row in calls:
                self.assertEqual(len(row), 24)
                self.assertEqual(row[10], "0")
            return [tuple(map(int, row[11:16])) for row in calls]

    def assert_support_partition(self, counts):
        for total, through, left, right, segmental in counts:
            with self.subTest(counts=(total, through, left, right, segmental)):
                self.assertTrue(all(count >= 0 for count in (total, through, left, right, segmental)))
                self.assertEqual(total, through + left + right + segmental)

    def test_exact_duplicate_alignment_counts_once(self):
        for tsd_enabled in (False, True):
            for sva in (False, True):
                with self.subTest(tsd_enabled=tsd_enabled, sva=sva):
                    counts = self.support_counts([SupportRow(0), SupportRow(0)], tsd_enabled, sva)
                    self.assert_support_partition(counts)
                    self.assertEqual(counts, [(1, 1, 0, 0, 0)])

    def test_coordinate_hsps_from_one_alignment_do_not_inflate_counts(self):
        rows = [SupportRow(0), SupportRow(0, 2, 299, 101, 399)]
        for tsd_enabled in (False, True):
            with self.subTest(tsd_enabled=tsd_enabled):
                counts = self.support_counts(rows, tsd_enabled)
                self.assert_support_partition(counts)
                # Existing clustering may retain separate coordinate variants.
                # Each emitted cluster must count the shared alignment once.
                self.assertTrue(all(count == (1, 1, 0, 0, 0) for count in counts))

    def test_two_independent_alignments_remain_two_supports(self):
        for tsd_enabled in (False, True):
            with self.subTest(tsd_enabled=tsd_enabled):
                counts = self.support_counts([SupportRow(0), SupportRow(1)], tsd_enabled)
                self.assert_support_partition(counts)
                self.assertEqual(counts, [(2, 2, 0, 0, 0)])

    def test_duplicate_seed_that_becomes_left_counts_once(self):
        left = SupportRow(0, 1, 100, 100, 200)
        for tsd_enabled in (False, True):
            with self.subTest(tsd_enabled=tsd_enabled):
                counts = self.support_counts([left, left, SupportRow(1)], tsd_enabled)
                self.assert_support_partition(counts)
                self.assertEqual(counts, [(2, 1, 1, 0, 0)])

    def test_same_alignment_cannot_count_as_go_through_and_left(self):
        rows = [SupportRow(0, 1, 100, 100, 200), SupportRow(0), SupportRow(1)]
        for tsd_enabled in (False, True):
            with self.subTest(tsd_enabled=tsd_enabled):
                counts = self.support_counts(rows, tsd_enabled)
                self.assert_support_partition(counts)
                self.assertIn((2, 1, 1, 0, 0), counts)
                self.assertTrue(all(count in {(2, 1, 1, 0, 0), (1, 1, 0, 0, 0)} for count in counts))

    def test_mixed_support_categories_stay_nonnegative(self):
        left = SupportRow(2, 1, 100, 100, 200)
        partial = SupportRow(3, 200, 300, 300, 400)
        rows = [SupportRow(0), SupportRow(0), SupportRow(1), left, left, partial, partial]
        for tsd_enabled in (False, True):
            with self.subTest(tsd_enabled=tsd_enabled):
                counts = self.support_counts(rows, tsd_enabled)
                self.assert_support_partition(counts)
                # Preserve the existing geometry: the terminal partial HSP in
                # this fixture is segmental rather than right-end support.
                self.assertIn((4, 2, 1, 0, 1), counts)
                self.assertTrue(all(count in {(4, 2, 1, 0, 1), (1, 1, 0, 0, 0)} for count in counts))

    def test_missing_tsd_id_does_not_create_go_through_support(self):
        for tsd_enabled in (False, True):
            with self.subTest(tsd_enabled=tsd_enabled):
                counts = self.support_counts([SupportRow(0)], tsd_enabled, keep_tsd_ids=False)
                self.assert_support_partition(counts)
                self.assertEqual(counts, [(1, 0, 0, 0, 1)])


if __name__ == "__main__":
    unittest.main()
