from CoolProp.CoolProp import PropsSI
import numpy as np

'''
Compressor

**Parameters**
    power: *float*
        Power of the compressor [W]
    eta_H: *float, optional*
        Enthalpic efficiency [0 to 1]
    eta_S: *float, optional*
        Entropic efficiency [0 to 1]
'''
class Compressor:
    def __init__(self, power=0, eta_H=1, eta_S=1):
        self.power = power
        self.eta_H = eta_H
        self.eta_S = eta_S
        
    def compress(self, fluid):
        if int(fluid.q) is not -1:
            h1 = PropsSI('H', 'T', fluid.t, 'Q', fluid.q, fluid.comp)
            s1 = PropsSI('S', 'T', fluid.t, 'Q', fluid.q, fluid.comp)
        else:
            h1 = PropsSI('H', 'T', fluid.t, 'P', fluid.p, fluid.comp)
            s1 = PropsSI('S', 'T', fluid.t, 'P', fluid.p, fluid.comp)
        
        h2 = h1 + (self.eta_H * self.power / fluid.flow)
        h2_isentropic = h1 + (self.eta_H * self.eta_S * self.power / fluid.flow)
        
        p2 = PropsSI('P', 'H', h2_isentropic, 'S', s1, fluid.comp)
        t2 = PropsSI('T', 'H', h2, 'P', p2, fluid.comp)
        q2 = PropsSI('Q', 'H', h2, 'P', p2, fluid.comp)
        
        fluid.set_state(new_temp=t2, new_pressure=p2, new_quality=q2)
        
        return fluid


'''
Joule-Thomson Valve

**Parameters**
    pressure_out: *float*

Currently assumes isenthalpic    

'''
class JT_Valve:
    def __init__(self,pressure_out=101325):
        self.pressure_out = pressure_out
    
    def expand(self, fluid):
        if fluid.q is not -1:
            h1 = PropsSI('H', 'T', fluid.t, 'Q', fluid.q, fluid.comp)
        else:
            h1 = PropsSI('H', 'T', fluid.t, 'P', fluid.p, fluid.comp)
            
        new_temp = PropsSI('T', 'H', h1, 'P', self.pressure_out, fluid.comp)
        new_quality = PropsSI('Q', 'H', h1, 'P', self.pressure_out, fluid.comp)
        
        fluid.set_state(new_temp=new_temp, new_pressure=self.pressure_out, new_quality=new_quality)
        
        return fluid



'''
Condenser
    Completely condenses a fluid at its current pressure

**Parameters**
    cond_press_drop: (float, optional)
        Pressure drop for the condensing fluid [Pa]
    cold_press_drop: (float, optional)
        Pressure drop for the cold fluid [Pa]

'''
class Condenser:
    def __init__(self, cond_press_drop=0, cold_press_drop=0):
        self.cond_press_drop = cond_press_drop
        self.cold_press_drop = cold_press_drop
    
    def exchange(self, cond_fluid, cold_fluid):
        
        assert cond_fluid.t > cold_fluid.t, \
            "Condensing fluid temeprature is less than the cold fluid temperature."
        
        # Required in case the quality is near 0 or 1
        if int(cond_fluid.q) is not -1:
            h_cond1 = PropsSI('H', 'P', cond_fluid.p, 'Q', cond_fluid.q, cond_fluid.comp)
        else:
            h_cond1 = PropsSI('H', 'P', cond_fluid.p, 'T', cond_fluid.t, cond_fluid.comp)
        
        # Drop quality to 0 for saturated liquid
        p_cond2 = cond_fluid.p - self.cond_press_drop
        h_cond2 = PropsSI('H', 'P', p_cond2, 'Q', 0, cond_fluid.comp)
        t_cond2 = PropsSI('T', 'P', p_cond2, 'Q', 0, cond_fluid.comp)
        cond_fluid.set_state(new_temp=t_cond2, new_pressure=p_cond2, new_quality=0)
        
        # Transfer the heat to the hot fluid
        q_transfer = (h_cond2 - h_cond1) * cond_fluid.flow
        
        h_cold1 = PropsSI('H', 'T', cold_fluid.t, 'P', cold_fluid.p, cold_fluid.comp)   
        h_cold2 = h_cold1 - q_transfer/cold_fluid.flow
        
        p_cold2 = cold_fluid.p - self.cold_press_drop
        t_cold2 = PropsSI('T', 'H', h_cold2, 'P', p_cold2, cold_fluid.comp)
        
        q_cold2 = PropsSI('Q', 'H', h_cold2, 'P', p_cold2, cold_fluid.comp)
        
        cold_fluid.set_state(new_temp=t_cold2, new_pressure=p_cold2, new_quality=q_cold2)
        
        assert t_cond2 > t_cold2, \
            "Impossible temperature crossover in Condenser."
        
        return cond_fluid, cold_fluid


'''
Evaporator
    Completely evaporates a fluid at its current pressure

**Parameters**
    evap_press_drop: (float, optional)
        Pressure drop for the evaporating fluid [Pa]
    hot_press_drop: (float, optional)
        Pressure drop for the hot fluid [Pa]
'''
class Evaporator:
    def __init__(self, evap_press_drop=0, hot_press_drop=0):
        self.evap_press_drop = evap_press_drop
        self.hot_press_drop = hot_press_drop
    
    def exchange(self, evap_fluid, hot_fluid):
        
        assert evap_fluid.t < hot_fluid.t, \
            "Evaporating fluid temeprature is greater than the hot fluid temperature."
        
        # Required in case the quality is near 0 or 1
        if evap_fluid.q is not -1:
            h_evap1 = PropsSI('H', 'P', evap_fluid.p, 'Q', evap_fluid.q, evap_fluid.comp)
        else:
            h_evap1 = PropsSI('H', 'P', evap_fluid.p, 'T', evap_fluid.t, evap_fluid.comp)
        
        # Increase quality to 1 for saturated vapor
        p_evap2 = evap_fluid.p - self.evap_press_drop
        h_evap2 = PropsSI('H', 'P', p_evap2, 'Q', 1, evap_fluid.comp)
        t_evap2 = PropsSI('T', 'P', p_evap2, 'Q', 1, evap_fluid.comp)
        evap_fluid.set_state(new_temp=t_evap2, new_pressure=p_evap2, new_quality=1)
        
        # Transfer the heat to the cold fluid
        q_transfer = (h_evap2 - h_evap1) * evap_fluid.flow
        
        h_hot1 = PropsSI('H', 'T', hot_fluid.t, 'P', hot_fluid.p, hot_fluid.comp)   
        h_hot2 = h_hot1 - q_transfer/hot_fluid.flow
        
        p_hot2 = hot_fluid.p - self.hot_press_drop
        t_hot2 = PropsSI('T', 'H', h_hot2, 'P', p_hot2, hot_fluid.comp)
        q_hot2 = PropsSI('Q', 'H', h_hot2, 'P', p_hot2, hot_fluid.comp)
        hot_fluid.set_state(new_temp=t_hot2, new_pressure=p_hot2, new_quality=q_hot2)
        
        assert t_evap2 < t_hot2, \
            "Impossible temperature crossover in Evaporator."
        
        return evap_fluid, hot_fluid


'''
Heat Exchanger
    Exchanges heat between two fluids

**Parameters**
    UA: (float)
        Overall heat transfer coefficient [J/K]
    hot_press_drop: (float, optional)
        Pressure drop for the hot fluid [Pa]
    cold_press_drop: (float, optional)
        Pressure drop for the cold fluid [Pa]
    hx_max_iter: (int, optional)
        Maximum number of iterations to solve for heat exchanger convergence
    hx_tol: (float, optional)
        Tolerance to stop heat exchanger convergence calculation [%]
'''
class Heat_Exchanger:
    def __init__(self, UA, hot_press_drop=0, cold_press_drop=0, hx_max_iter=50, hx_tol=0.01):
        self.UA = UA
        self.hot_press_drop = hot_press_drop
        self.cold_press_drop = cold_press_drop
        self.hx_max_iter = hx_max_iter
        self.hx_tol = hx_tol

    def exchange(self, hot_fluid, cold_fluid):
        
        if hot_fluid.q is not -1:
            h_hot1 = PropsSI('H', 'P', hot_fluid.p, 'Q', hot_fluid.q, hot_fluid.comp)
        else:
            h_hot1 = PropsSI('H', 'P', hot_fluid.p, 'T', hot_fluid.t, hot_fluid.comp)
        
        if cold_fluid.q is not -1:
            h_cold1 = PropsSI('H', 'P', cold_fluid.p, 'Q', cold_fluid.q, cold_fluid.comp)
        else:
            h_cold1 = PropsSI('H', 'P', cold_fluid.p, 'T', cold_fluid.t, cold_fluid.comp)
        
        # Put initial temperatures in consistently named variables
        t_hot1 = hot_fluid.t
        t_cold1 = cold_fluid.t
        
        # Adjust final pressures
        p_hot2 = hot_fluid.p - self.hot_press_drop
        p_cold2 = cold_fluid.p - self.cold_press_drop
        
        # Iteratively calculate outlet temperatures
        q_guess = self.UA * (t_hot1 - t_cold1)
        for i in range(self.hx_max_iter):
            h_hot2 = h_hot1 - (q_guess / hot_fluid.flow)
            h_cold2 = h_cold1 + (q_guess / cold_fluid.flow)
            
            t_hot2 = PropsSI('T', 'H', h_hot2, 'P', p_hot2, hot_fluid.comp)
            t_cold2 = PropsSI('T', 'H', h_cold2, 'P', p_cold2, cold_fluid.comp)
            
            dT1 = t_hot1 - t_cold2
            dT2 = t_hot2 - t_cold1
            
            dTlm = (dT2 - dT1) / np.log(dT2/dT1)
            
            
            print(dT1, dT2, dTlm, q_guess)
            
            
            q_actual = self.UA * dTlm
            pcent_diff = 2*abs(q_guess - q_actual) / (q_guess + q_actual)
            
            if pcent_diff < self.hx_tol:
                break
            else:
                q_guess = q_actual
        else:
            print("Heat exchanger did not converge within the specified tolerance")
        
        
        
        q_hot2 = PropsSI('Q', 'H', h_hot2, 'P', p_hot2, hot_fluid.comp)
        q_cold2 = PropsSI('Q', 'H', h_cold2, 'P', p_cold2, cold_fluid.comp)
        
        hot_fluid.set_state(new_temp=t_hot2, new_pressure=p_hot2, new_quality=q_hot2)
        cold_fluid.set_state(new_temp=t_cold2, new_pressure=p_cold2, new_quality=q_cold2)
        
        return hot_fluid, cold_fluid


'''
Fluid object

**Parameters**
    comp: *string*
        Name of the first component
    flow: *float*
        Fluid flowrate [kg/s]
    temperature: *float, optional*
        Fluid temperature [K]
    pressure: *float, optional*
        Fluid pressure [Pa]
    quality: *float, optional*
        Vapor quality [0 to 1]
'''
class Pure_Fluid:
    def __init__(self, comp, flow, temperature=None, pressure=None, quality=None):
        self.comp = comp
        self.flow = flow
        
        assert (temperature,pressure,quality).count(None) == 1, \
            "Must specify exactly two properties for a pure fluid."
        
        if temperature is not None and quality is not None:
            self.t = temperature
            self.q = quality
            self.p = PropsSI('P', 'T', temperature, 'Q', quality, comp)
        elif pressure is not None and quality is not None:
            self.p = pressure
            self.q = quality
            self.t = PropsSI('T', 'P', pressure, 'Q', quality, comp)
        elif temperature is not None and pressure is not None:
            self.t = temperature
            self.p = pressure
            self.q = -1

    
    def set_state(self, new_temp=None, new_pressure=None, new_quality=None):
        if new_temp is not None:
            self.t = new_temp
        if new_pressure is not None:
            self.p = new_pressure
        if new_quality is not None:
            self.q = new_quality


# 8.3 gpm is 0.52 kg/s 294K

cond_press = 9 * 101325

f = Pure_Fluid(comp='R410A', flow=0.07, pressure = cond_press, quality = 1)
w = Pure_Fluid(comp='Water', flow=0.52, temperature = 300, pressure = 2 * 101325)
a = Pure_Fluid(comp='air', flow=0.58, temperature = 350, pressure = 1.1 * 101325)

'''
comp = Compressor(power=2250, eta_H=0.98, eta_S=0.70)
cond = Condenser(cond_press_drop=5000, cold_press_drop=5000)
valve = JT_Valve(pressure_out = cond_press)
evap = Evaporator(evap_press_drop=50000, hot_press_drop=10000)

print(int(f.t), round(f.p/101325,2), round(f.q,2))
f = comp.compress(f)
print(int(f.t), round(f.p/101325,2), round(f.q,2))
f,w = cond.exchange(f,w)
print(int(f.t), round(f.p/101325,2), round(f.q,2))
f = valve.expand(f)
print(int(f.t), round(f.p/101325,2), round(f.q,2))
f,a = evap.exchange(f,a)
print(int(f.t), round(f.p/101325,2), round(f.q,2))

print('w', int(w.t), round(w.p/101325,2), round(w.q,2))
print('a', int(a.t), round(a.p/101325,2), round(a.q,2))
'''


hx = Heat_Exchanger(UA=570)

print('air', float(a.t), 'water', float(w.t))
a,w = hx.exchange(a,w)
print('air', round(a.t,1), 'water', round(w.t,1))





## More Space

