import numpy as np

class Counters:
    def __init__(self):
        self._nWithinAcceptance = 0
        self._nWithinTPC = 0
        self._twoTracks = 0
        self._twoTruthTracks = 0
        self._allTracks = 0
        self._nMatchingTracks = 0
        self._nMatchingTracksReffited = 0
        self._nMatchingVertices = 0
        self._matchingVerticesTPC = 0
        self._twoRecoLeptons = 0
        self._onlyMarlinTrkTracks = 0
        self._onlyRefittedTracks = 0
        self._nPFOs = 0
        self._nCorrPID = 0
        self._nEvSecVertices = 0
        self._nEvSecVerticesTPC = 0
        self._nSecVertices = 0
        self._nSecVerticesTPC = 0
        self._nFakeVertices = 0
        self._nFakeVerticesTPC = 0
        self._nPhotonConversions = 0
        self._nPhotonConversionsTPC = 0
        self._secondaryIntCut = 0
        self._eeMassCut = 0
        self._pipiMassCut = 0
        self._lambdaMassCut = 0
        self._pointingIPCut = 0
        self._weightedEvents = np.zeros(9)
        self._weightedEventsInTPC = np.zeros(9)
        self._matchingWeightedEventsTPC = np.zeros(9)
        self._weightedEventsErr = np.zeros(9)
        self._matchingWeightedEventsTPCErr = np.zeros(9)

    def increment(self, variable, condition, value):
        if condition:
            setattr(self, variable, getattr(self, variable) + value)