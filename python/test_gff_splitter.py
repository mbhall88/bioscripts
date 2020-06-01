import pytest
from intervaltree import IntervalTree, Interval

from gff_splitter import merge_short_intervals, infer_interval_type, NoDataError


class TestMergeShortIntervals:
    def test_emptyTree_returnsEmptyTree(self):
        tree = IntervalTree()

        actual = merge_short_intervals(tree)

        assert actual.is_empty()

    def test_oneInterval_returnsSameTree(self):
        tree = IntervalTree([Interval(1, 3)])
        min_len = 1

        actual = merge_short_intervals(tree, min_len=min_len)
        expected = tree

        assert actual == expected

    def test_oneIntervalBelowMinLen_returnIntervalEndPlusOne(self):
        tree = IntervalTree([Interval(1, 3)])
        min_len = 9

        actual = merge_short_intervals(tree, min_len=min_len)
        expected = IntervalTree([Interval(1, 4)])

        assert actual == expected

    def test_twoIntervalsAboveMinLen_returnSame(self):
        tree = IntervalTree([Interval(1, 3), Interval(3, 7)])
        min_len = 1

        actual = merge_short_intervals(tree, min_len=min_len)
        expected = tree

        assert actual == expected

    def test_twoIntervalsOneBelowMinLen_returnOneInterval(self):
        tree = IntervalTree([Interval(1, 3, "foo"), Interval(3, 7, "bar")])
        min_len = 3

        actual = merge_short_intervals(tree, min_len=min_len)
        expected = IntervalTree([Interval(1, 7, data="foo+bar")])

        assert actual == expected

    def test_threeIntervalsAllBelowMinLen_returnOneInterval(self):
        tree = IntervalTree(
            [Interval(1, 3, "foo"), Interval(3, 7, "bar"), Interval(7, 11, "baz")]
        )
        min_len = 9

        actual = merge_short_intervals(tree, min_len=min_len)
        expected = IntervalTree([Interval(1, 11, data="foo+bar+baz")])

        assert actual == expected

    def test_threeIntervalsTwoBelowMinLen_returnTwoIntervals(self):
        tree = IntervalTree(
            [Interval(1, 3, "foo"), Interval(3, 7, "bar"), Interval(7, 19, "baz")]
        )
        min_len = 6

        actual = merge_short_intervals(tree, min_len=min_len)
        expected = IntervalTree(
            [Interval(1, 7, data="foo+bar"), Interval(7, 19, "baz")]
        )

        assert actual == expected


class TestInferIntervalType:
    def test_intervalWithNoData_raisesError(self):
        interval = Interval(1, 2)
        with pytest.raises(NoDataError):
            infer_interval_type(interval)

    def test_intervalIsFeature_returnsFeature(self):
        interval = Interval(1, 2, data="rpoB")

        actual = infer_interval_type(interval)
        expected = "feature"

        assert actual == expected

    def test_intervalIsMergedFeature_returnsMergedFeature(self):
        interval = Interval(1, 20, data="rpoB+rpoC")

        actual = infer_interval_type(interval)
        expected = "merged_features"

        assert actual == expected

    def test_intervalIsIGR_returnsIGR(self):
        interval = Interval(1, 20, data="IGR:1-20")

        actual = infer_interval_type(interval)
        expected = "igr"

        assert actual == expected

    def test_intervalIsMergedIGRAndFeature_returnsMergedIGRAndFeature(self):
        interval = Interval(1, 20, data="rpoB+rpoC+IGR:19-20")

        actual = infer_interval_type(interval)
        expected = "merged_feature_and_igr"

        assert actual == expected

    def test_intervalIsMergedIGRs_returnsMergedIGRs(self):
        interval = Interval(1, 20, data="IGR:1-19+IGR:19-20")

        actual = infer_interval_type(interval)
        expected = "merged_igrs"

        assert actual == expected
