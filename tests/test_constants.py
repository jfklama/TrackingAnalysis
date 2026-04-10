import constants


def test_tpc_geometry_positive():
    assert constants.TPC_R_MAX > 0
    assert constants.TPC_Z_MAX > 0


def test_tpc_z_longer_than_r():
    # ILD is a barrel detector: half-length exceeds radius
    assert constants.TPC_Z_MAX > constants.TPC_R_MAX


def test_particle_masses_positive():
    assert constants.PROTON_MASS > 0
    assert constants.PION_MASS > 0
    assert constants.LAMBDA_MASS > 0
    assert constants.K0S_MASS > 0


def test_mass_ordering():
    # PDG values: Lambda(1115) > proton(938) > pion(140) MeV
    assert constants.LAMBDA_MASS > constants.PROTON_MASS > constants.PION_MASS


def test_k0s_lighter_than_proton():
    assert constants.K0S_MASS < constants.PROTON_MASS


def test_electron_mass_zero():
    # Ultra-relativistic approximation used throughout the analysis
    assert constants.ELECTRON_MASS == 0.0


def test_nbins_positive():
    assert constants.NBINS > 0
