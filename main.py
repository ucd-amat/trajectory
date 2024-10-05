from APCreader import parse_APC_data
import numpy as np
import matplotlib.pyplot as plt

# Usage
file_name = 'data/PER3_16x8.dat'  # Replace with your actual filename
prop_performance = parse_APC_data(file_name)[0]
# print(prop_performance(2,80))
thrust_lbf = 2
airspeed_fps = 80
# Power_W = prop_performance(thrust_lbf,airspeed_fps)["PWR_W"]
Weight = 10 #lbs
S = 6 #ft^2 
b = 6 #ft
AR = 6
rho = 0.002387 #slugs/ft^3


def dragCoeff_estimate(CL=0.5,CD0=0.052,AR=6,e=0.7):
    CD = CD0 + CL**2/(np.pi*e*AR)
    return CD

def cruise(distance,airspeed_fps):
    qS = ((1/2)*rho*airspeed_fps**2*S)
    CL = Weight/qS
    print('Cruise CL', CL)
    CD = dragCoeff_estimate(CL)
    D = CD*qS
    thrust_lbf = D # assuming drag = thrust
    Power_W = prop_performance(thrust_lbf,airspeed_fps)["PWR_W"]
    Energy_J = Power_W*(distance/airspeed_fps) #energy used by propeller estimate
    
    #convert joules to Watts*hours (1 J = 1/3600 Wh)
    eta_Motor = 0.8 # motor efficiency
    Watt_hours = (Energy_J/eta_Motor)/3600 # convert to Watt*Hours
    return Watt_hours

def half_turn(n,airspeed_fps):
    qS = ((1/2)*rho*airspeed_fps**2*S)
    CL = (Weight*n)/qS
    CD = dragCoeff_estimate(CL)
    D = CD*qS
    phi = np.arccos(1/n) # get bank angle from load factor
    r = airspeed_fps**2/(32.2*np.tan(phi)) # radius of turn
    distance = np.pi*r # distance for half turn
    thrust_lbf = D # assuming drag = thrust
    Power_W = prop_performance(thrust_lbf,airspeed_fps)["PWR_W"]
    Energy_J = Power_W*(distance/airspeed_fps) #energy used by propeller estimate
    
    #convert joules to Watts*hours (1 J = 1/3600 Wh)
    eta_Motor = 0.8 # motor efficiency
    Watt_hours = (Energy_J/eta_Motor)/3600 # convert to Watt*Hours
    phi_deg = phi*180/np.pi
    return Watt_hours, CL, phi_deg

# Lap Simulation
distance = 500*4 #ft
airspeed_fps = 80
n = 1.4 # load factor

# total energy spent in Watt Hours at constant altitude for 1 lap
WattHours_tot = cruise(distance,airspeed_fps) + 4*half_turn(n,airspeed_fps)[0]

print('Turn CL',half_turn(n,airspeed_fps)[1])

print('Total WattHours Spent',WattHours_tot)

def takeoff(V_rot,theta_climb_deg,rotation_rate_degps,gtow,RPM = 10000):
    x = np.zeros(4) #initialize x
    rotating = False
    t = 0
    dt = 0.01
    theta = 0 #radians
    theta0 = 2*(np.pi/180) #radius
    rotation_rate = rotation_rate_degps*(np.pi/180)
    prop_performance = parse_APC_data(file_name,x = "RPM", y = "V_fps")[0]
    theta_climb = theta_climb_deg*(np.pi/180) # deg to rads

    history = dict(
        t = [],
        x = [],
        y = [],
        u = [],
        v = [],
        theta = [],
        alpha = [],
        CL = [],
        CD = [],
        Power = [],

    )

    # while x[2] < 35 and x[0] < 1000: #ft
    while x[0] < 1000: #ft
        t += dt
        if x[1] > V_rot and x[2] == 0 and not rotating:
            rotating = True
            t_rotation = t

        if rotating:
            # rotate until you hit the desired climb angle
            theta = min((t-t_rotation)*rotation_rate,theta_climb)

        normal_airspeed_fps = ((x[1]*np.cos(theta))**2 + (x[3]*np.sin(theta))**2)**(0.5)
        airspeed_fps = ((x[1])**2 + (x[3])**2)**(0.5)

        prop_data = prop_performance(RPM,normal_airspeed_fps)
        Thrust = prop_data["Thrust_lbf"]
        Power = prop_data["PWR_W"]
        qS = ((1/2)*rho*airspeed_fps**2*S)
        alpha = theta - np.nan_to_num(np.arctan(x[3]/x[1]))
        CL = 2*np.pi*(alpha + theta0) #assume linear lift slope
        CD = dragCoeff_estimate(CL)
        D = CD*qS
        L = CL*qS

        dx = np.array([
            x[1], #x_dot
            (Thrust*np.cos(theta)-D*np.cos(theta)-L*np.sin(theta))/(gtow/32.2),
            x[3], #y_dot
            (Thrust*np.sin(theta)-D*np.sin(theta)+L*np.cos(theta))/(gtow/32.2) - 32.2,
        ])

        if dx[3] < 0 and x[2] <= 0: #y_2dot < 0
            dx[2] = 0 #y_dot
            dx[3] = 0 #y_2dot

        x += dx*dt

        # print("Time (s)", t)
        # print("Theta",theta*(180/np.pi))
        # print("Alpha",alpha*(180/np.pi))
        # print("dx",dx)
        # print("x",x)
        # print("CL",CL)

        history["t"].append(t)
        history["x"].append(x[0])
        history["y"].append(x[2])
        history["u"].append(x[1])
        history["v"].append(x[3])
        history["theta"].append(theta)
        history["alpha"].append(alpha)
        history["CL"].append(CL)
        history["CD"].append(CD)
        history["Power"].append(Power)

    return history
        

history = takeoff(50,15,10,15,8000)

for var in ["x","y","u","v","theta","alpha","CL","CD","Power"]:
    if var == "theta" or var == "alpha":
        plt.plot(history["t"],np.array(history[var])*(180/np.pi),label=var)
    else:
        plt.plot(history["t"],history[var],label=var)

plt.xlabel("t")
plt.legend()
plt.grid()

plt.figure()

for var in ["t","y","u","v","theta","alpha","CL","CD","Power"]:
    if var == "theta" or var == "alpha":
        plt.plot(history["x"],np.array(history[var])*(180/np.pi),label=var)
    else:
        plt.plot(history["x"],history[var],label=var)

plt.xlabel("x")
plt.legend()
plt.grid()

print("Takeoff Energy (Wh)", np.trapz(history["Power"],x=history["t"])/3600)

plt.show()


    





    
