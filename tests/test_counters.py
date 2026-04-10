import numpy as np
import pytest

from counters import Counters


INTEGER_COUNTERS = [
    '_nWithinAcceptance', '_nWithinTPC', '_twoTracks', '_twoTruthTracks',
    '_allTracks', '_nMatchingTracks', '_nMatchingTracksReffited',
    '_nMatchingVertices', '_matchingVerticesTPC', '_twoRecoLeptons',
    '_onlyMarlinTrkTracks', '_onlyRefittedTracks', '_nPFOs', '_nCorrPID',
    '_nEvSecVertices', '_nEvSecVerticesTPC', '_nSecVertices',
    '_nSecVerticesTPC', '_nFakeVertices', '_nFakeVerticesTPC',
    '_nPhotonConversions', '_nPhotonConversionsTPC', '_secondaryIntCut',
    '_eeMassCut', '_pipiMassCut', '_lambdaMassCut', '_pointingIPCut',
]

ARRAY_COUNTERS = [
    '_weightedEvents', '_weightedEventsInAcceptance', '_weightedEventsInTPC',
    '_matchingWeightedEventsAcceptance', '_matchingWeightedEventsTPC',
    '_weightedEventsErr', '_matchingWeightedEventsAcceptErr',
    '_matchingWeightedEventsTPCErr',
]


def test_init_integer_counters_zero():
    c = Counters()
    for attr in INTEGER_COUNTERS:
        assert getattr(c, attr) == 0, f"{attr} should be 0 on init"


def test_init_numpy_arrays():
    c = Counters()
    for attr in ARRAY_COUNTERS:
        arr = getattr(c, attr)
        assert isinstance(arr, np.ndarray), f"{attr} should be ndarray"
        assert len(arr) == 9, f"{attr} should have length 9"
        assert np.all(arr == 0), f"{attr} should be all zeros"


def test_increment_when_true():
    c = Counters()
    c.increment('_nMatchingTracks', True, 1)
    assert c._nMatchingTracks == 1


def test_no_increment_when_false():
    c = Counters()
    c.increment('_nMatchingTracks', False, 1)
    assert c._nMatchingTracks == 0


def test_increment_by_custom_value():
    c = Counters()
    c.increment('_nSecVertices', True, 3)
    assert c._nSecVertices == 3


def test_multiple_increments_accumulate():
    c = Counters()
    for _ in range(5):
        c.increment('_twoTracks', True, 1)
    assert c._twoTracks == 5


def test_increment_mixed_conditions():
    c = Counters()
    c.increment('_nMatchingVertices', True, 1)
    c.increment('_nMatchingVertices', False, 1)
    c.increment('_nMatchingVertices', True, 1)
    assert c._nMatchingVertices == 2
