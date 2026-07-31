from euv_playground.resist import shot_noise


def test_reference_shot_noise():
    _, photons, noise = shot_noise(30, 10)
    assert 15 <= photons <= 25
    assert 0.18 <= noise <= 0.27
