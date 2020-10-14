# -*- coding: utf-8 -*-

import thermodynamic_toolbox as thermo

mix = {"N2": 0.78084, "O2": 0.20942, "Ar": 0.00934, "CO2": 0.0004}

def equal_with_margin(val1, val2, margin=0.00001):
    print("Within " + str(abs(val1 - val2) / val2 * 1000000) + " ppm, acceptaple " + str(margin * 1000000) + " ppm")
    if abs(val1 - val2) / val2 <= margin:
        return True
    return False

def test_equal_with_margin_true():
    assert equal_with_margin(2, 2)
    
def test_equal_with_margin_false():
    assert equal_with_margin(3, 3.0001) == False

def test_isentrop_press_temp():
    assert equal_with_margin(thermo.isentrop_press_temp(1.3, 100000, 300, 330), 151136.1317039667157804610127936704129735)

def test_isentrop_press_vol():
    assert equal_with_margin(thermo.isentrop_press_vol(1.3, 100000, 1, 1.3), 71100.6616149676017537443064557694)

def test_isentrop_temp_press():
    assert equal_with_margin(thermo.isentrop_temp_press(1.3, 300, 100000, 150000), 329.425868143664876303540348645287182)

def test_isentrop_temp_vol():
    assert equal_with_margin(thermo.isentrop_temp_vol(1.3, 300, 1, 0.9), 309.6338992284570104848606920025972024101842)
    
def test_isentrop_vol_temp():
    assert equal_with_margin(thermo.isentrop_vol_temp(1.3, 1, 300, 330), 0.72782066577871095717418012108415351238)

def test_isentrop_vol_press():
    assert equal_with_margin(thermo.isentrop_vol_press(1.3, 1, 100000, 150000), 0.73205748476369972511897855254508262670640)
    
def test_mix_to_CP_string():
    assert thermo.mix_to_CP_string(mix) == "N2[0.78084]&O2[0.20942]&Ar[0.00934]&CO2[0.0004]"
    
def test_mass_frac():
    assert thermo.mass_frac(mix) == {'N2': 0.7551640362241803, 'O2': 0.23134708444921825, 'Ar': 0.012881134081446947, 'CO2': 0.0006077452451545295}
    
def test_molar_mixer():
    assert False
    
def test_mass_mixer():
    assert False
    
def test_get_kappa():
    assert False
    
def test_get_rho():
    assert False
    
def test_get_speed_of_sound():
    assert False