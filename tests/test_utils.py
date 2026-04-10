import math
from unittest.mock import MagicMock

import pytest

import utils


FCT = 2.99792458E-4
B_FIELD = 3.5  # hardcoded in utils.py


class TestGeometry:

    def test_R_pythagorean_triple(self):
        assert utils.R(3, 4) == pytest.approx(5.0)

    def test_R_zero(self):
        assert utils.R(0, 0) == 0.0

    def test_R_negative_inputs(self):
        assert utils.R(-3, -4) == pytest.approx(5.0)

    def test_R3_pythagorean(self):
        # 1^2 + 2^2 + 2^2 = 9 -> sqrt = 3
        assert utils.R3(1, 2, 2) == pytest.approx(3.0)

    def test_R3_zero(self):
        assert utils.R3(0, 0, 0) == 0.0

    def test_spacialDistance_known(self):
        v1 = [0, 0, 0]
        v2 = [3, 4, 0]
        assert utils.spacialDistance(v1, v2) == pytest.approx(5.0)

    def test_spacialDistance_symmetric(self):
        v1 = [1, 2, 3]
        v2 = [4, 6, 3]
        assert utils.spacialDistance(v1, v2) == pytest.approx(utils.spacialDistance(v2, v1))

    def test_spacialDistance_zero(self):
        v = [1.5, -2.3, 7.0]
        assert utils.spacialDistance(v, v) == pytest.approx(0.0)


class TestPhiDistance:

    def test_same_angle(self):
        assert utils.phiDistance(1.0, 1.0) == pytest.approx(0.0)

    def test_quarter_pi(self):
        assert utils.phiDistance(0, math.pi / 4) == pytest.approx(math.pi / 4)

    def test_half_pi(self):
        assert utils.phiDistance(0, math.pi / 2) == pytest.approx(math.pi / 2)

    def test_wrap_around(self):
        # |0 - 1.5π| = 1.5π > π  →  2π - 1.5π = 0.5π
        assert utils.phiDistance(0, 1.5 * math.pi) == pytest.approx(0.5 * math.pi)

    def test_symmetric(self):
        assert utils.phiDistance(0.3, 1.2) == pytest.approx(utils.phiDistance(1.2, 0.3))

    def test_result_in_range(self):
        for a, b in [(0, 0.1), (0, math.pi), (1.0, 5.0), (0, 2 * math.pi - 0.01)]:
            result = utils.phiDistance(a, b)
            assert 0.0 <= result <= math.pi


class TestArcLength:

    def test_negative_charge_forward(self):
        # q < 0: dPhi = phi2 - phi1 = pi/2 - 0 = pi/2 (positive, no wrap)
        assert utils.getArcLength(0, math.pi / 2, -1) == pytest.approx(math.pi / 2)

    def test_positive_charge_needs_wrap(self):
        # q > 0: dPhi = 0 - pi/2 = -pi/2 < 0  =>  += 2pi  => 3pi/2
        assert utils.getArcLength(0, math.pi / 2, 1) == pytest.approx(1.5 * math.pi)

    def test_zero_arc(self):
        assert utils.getArcLength(1.0, 1.0, -1) == pytest.approx(0.0)
        assert utils.getArcLength(1.0, 1.0, 1) == pytest.approx(0.0)

    def test_full_circle(self):
        # q < 0: dPhi = 2pi - 0 = 2pi
        assert utils.getArcLength(0, 2 * math.pi, -1) == pytest.approx(2 * math.pi)


class TestAngularDistance:

    def test_same_direction_zero(self):
        assert utils.angularDistance(0.5, 0.5, 1.0, 1.0) == pytest.approx(0.0)

    def test_phi_only_difference(self):
        # theta equal, sin term cancels, dist = |dPhi|
        theta = math.pi / 4
        dphi = math.pi / 4
        expected = dphi  # sin term is zero when theta diff is zero
        assert utils.angularDistance(theta, theta, 0.0, -dphi) == pytest.approx(expected)

    def test_phi_wraps_above_pi(self):
        # dPhi = pi + 0.1 > pi  →  dPhi -= 2pi  => -(pi - 0.1)
        theta = math.pi / 2
        d = utils.angularDistance(theta, theta, math.pi + 0.1, 0.0)
        # dPhi = pi+0.1 - 0 > pi → dPhi = pi+0.1 - 2pi = -(pi-0.1)
        assert d == pytest.approx(math.pi - 0.1)


class TestRotateAngle:

    def test_no_wrap(self):
        assert utils.rotateAngle(math.pi / 4, math.pi / 4) == pytest.approx(math.pi / 2)

    def test_wrap_negative(self):
        # pi/4 + (-pi) = -3pi/4  ->  += 2pi  => 5pi/4
        result = utils.rotateAngle(math.pi / 4, -math.pi)
        assert result == pytest.approx(5 * math.pi / 4)

    def test_zero_rotation(self):
        angle = 1.23
        assert utils.rotateAngle(angle, 0) == pytest.approx(angle)


class TestTrackMomentum:

    def _make_track(self, phi, tanLambda, omega):
        track = MagicMock()
        track.getPhi.return_value = phi
        track.getTanLambda.return_value = tanLambda
        track.getOmega.return_value = omega
        return track

    def test_forward_track(self):
        # phi=0, tanLambda=0, omega=1  =>  pxy = FCT*B*1, px=pxy, py=0, pz=0
        track = self._make_track(0.0, 0.0, 1.0)
        mom = utils.getTrackMomentum(track)
        pxy = FCT * B_FIELD * 1.0
        assert mom[0] == pytest.approx(pxy)
        assert mom[1] == pytest.approx(0.0, abs=1e-12)
        assert mom[2] == pytest.approx(0.0, abs=1e-12)

    def test_transverse_momentum_from_omega(self):
        # radius = 1/|omega|, so larger omega => smaller radius => smaller pxy
        track_high = self._make_track(0.0, 0.0, 2.0)
        track_low = self._make_track(0.0, 0.0, 0.5)
        mom_high = utils.getTrackMomentum(track_high)
        mom_low = utils.getTrackMomentum(track_low)
        assert abs(mom_high[0]) < abs(mom_low[0])

    def test_zero_omega_gives_zero_momentum(self):
        track = self._make_track(0.0, 0.0, 0.0)
        mom = utils.getTrackMomentum(track)
        assert mom == [0.0, 0.0, 0.0]

    def test_longitudinal_momentum_from_tanLambda(self):
        track = self._make_track(0.0, 1.0, 1.0)
        mom = utils.getTrackMomentum(track)
        pxy = FCT * B_FIELD
        assert mom[2] == pytest.approx(1.0 * pxy)


class TestTrackCharge:

    def _make_track(self, omega):
        track = MagicMock()
        track.getOmega.return_value = omega
        return track

    def test_positive_omega_positive_charge(self):
        assert utils.getTrackCharge(self._make_track(0.5)) == pytest.approx(1.0)

    def test_negative_omega_negative_charge(self):
        assert utils.getTrackCharge(self._make_track(-0.5)) == pytest.approx(-1.0)

    def test_large_omega_sign_only(self):
        assert utils.getTrackCharge(self._make_track(100.0)) == pytest.approx(1.0)
        assert utils.getTrackCharge(self._make_track(-100.0)) == pytest.approx(-1.0)
