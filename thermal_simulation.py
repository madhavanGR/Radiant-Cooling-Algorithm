
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --- Wall and Material Properties ---
Lx, Ly = 0.2, 2.5
k, rho, cp = 0.8, 2000, 1000
alpha = k / (rho * cp)

# --- Grid Setup ---
dx = dy = 0.01
nx, ny = int(Lx / dx) + 1, int(Ly / dy) + 1
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
dt = 60  # seconds

# --- Time Setup ---
total_hours = 24
time_steps = int((total_hours * 3600) / dt)

# --- Initial Conditions ---
T_initial, T_indoor = 20.0, 22.0
T = np.full((ny, nx), T_initial)

# --- Chilled Water Loop ---
T_water_in = 16.0
cp_water = 4186
flow_rate = 0.1
U_pipe = 200
pipe_spacing = 0.2
pipe_depth_index = nx - 3
pipe_columns = np.arange(0, ny, int(pipe_spacing / dy))

# --- Indoor Air Properties ---
room_width, room_height, room_depth = 5, Ly, 4
air_density, air_cp = 1.225, 1005
air_volume = room_width * room_height * room_depth
air_mass = air_density * air_volume

# --- HVAC Control ---
T_set_high, T_set_low = 24.0, 22.0
cooling_active = False

# --- Load Weather Data ---
df = pd.read_csv('bengaluru_hourly_temps.csv')
hourly_temps = df['Temperature'].values
hourly_humidity = df['Humidity'].values
external_temperature = np.interp(np.linspace(0, 23.99, time_steps), np.arange(24), hourly_temps)
external_humidity = np.interp(np.linspace(0, 23.99, time_steps), np.arange(24), hourly_humidity)

# --- Result Storage ---
T_indoor_list, T_water_out_list, cooling_status_list = [], [], []
pmv_list, ppd_list = [], []

# --- PMV/PPD Function ---
def calc_pmv(tdb, tr, vr, rh, met=1.1, clo=0.7):
    pa = rh * 10 * np.exp(16.6536 - 4030.183 / (tdb + 235))
    icl = 0.155 * clo
    m = met * 58.15
    w = 0
    mw = m - w
    fcl = 1.05 + 0.1 * clo if clo > 0.5 else 1 + 0.2 * clo
    hcf = 12.1 * np.sqrt(vr)
    taa = tdb + 273
    tra = tr + 273
    tcla = taa + (35.5 - tdb) / (3.5 * (6.45 * icl + 1.0))
    for _ in range(100):
        tcl_old = tcla
        hc = max(2.38 * abs(100.0 * (tcla - taa))**0.25, hcf)
        tcla = (35.7 - 0.028 * mw + (fcl * hc * taa + fcl * 3.96e-8 * (tra**4)) * icl) /                (1 + fcl * (hc + 3.96e-8 * ((tra**4 + taa**4))) * icl)
        if abs(tcla - tcl_old) < 0.001:
            break
    tcl = tcla - 273
    hl1 = 3.05 * 0.001 * (5733 - (6.99 * mw) - pa)
    hl2 = 0.42 * (mw - 58.15) if mw > 58.15 else 0
    hl3 = 1.7 * 0.00001 * m * (5867 - pa)
    hl4 = 0.0014 * m * (34 - tdb)
    hl5 = fcl * hc * (tcl - tdb)
    hl6 = fcl * 3.96e-8 * ((tcl + 273)**4 - (tra)**4)
    ts = mw - (hl1 + hl2 + hl3 + hl4 + hl5 + hl6)
    pmv = (0.303 * np.exp(-0.036 * m) + 0.028) * ts
    ppd = 100 - 95 * np.exp(-0.03353 * pmv**4 - 0.2179 * pmv**2)
    return pmv, ppd

# --- Heat Diffusion Function ---
def update_temperature(T, alpha, dx, dy, dt):
    T_new = T.copy()
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            d2Tdx2 = (T[j, i-1] - 2*T[j, i] + T[j, i+1]) / dx**2
            d2Tdy2 = (T[j-1, i] - 2*T[j, i] + T[j+1, i]) / dy**2
            T_new[j, i] += alpha * dt * (d2Tdx2 + d2Tdy2)
    return T_new

# --- Time Loop ---
for t in range(time_steps):
    T = update_temperature(T, alpha, dx, dy, dt)
    T[:, 0] = external_temperature[t]
    T[:, -1] = T_indoor

    if T_indoor > T_set_high:
        cooling_active = True
    elif T_indoor < T_set_low:
        cooling_active = False

    # --- Chilled Water Cooling ---
    if cooling_active:
        T_water = T_water_in
        total_heat_removed = 0
        for row in pipe_columns:
            wall_temp = T[row, pipe_depth_index]
            delta_T = wall_temp - T_water
            q_conv = U_pipe * delta_T
            heat_exchange = q_conv * dx * dy
            T[row, pipe_depth_index] -= (heat_exchange * dt) / (rho * cp * dx * dy)
            T_water += (heat_exchange * dt) / (flow_rate * cp_water)
            total_heat_removed += heat_exchange
    else:
        total_heat_removed = 0
        T_water = T_water_out_list[-1] if t > 0 else T_water_in

    inner_wall = T[:, -2]
    q_in = k * (inner_wall - T_indoor) / dx
    total_q_in = np.sum(q_in) * dy
    net_q = total_q_in - total_heat_removed
    dT_indoor = (net_q * dt) / (air_mass * air_cp)
    T_indoor += dT_indoor
    T[:, -1] = T_indoor

    # --- PMV/PPD Calculation ---
    Tr = T_indoor
    Va = 0.1
    Rh = external_humidity[t]
    pmv, ppd = calc_pmv(tdb=T_indoor, tr=Tr, vr=Va, rh=Rh)

    # --- Store ---
    T_indoor_list.append(T_indoor)
    T_water_out_list.append(T_water)
    cooling_status_list.append(cooling_active)
    pmv_list.append(pmv)
    ppd_list.append(ppd)

# --- Plotting ---
time_array = np.arange(time_steps) * dt / 3600

plt.figure(figsize=(10, 5))
plt.plot(time_array, T_indoor_list, label='Indoor Temp')
plt.axhline(T_set_high, color='r', linestyle='--', label='Setpoint High')
plt.axhline(T_set_low, color='b', linestyle='--', label='Setpoint Low')
plt.xlabel('Time (hours)')
plt.ylabel('Temperature (°C)')
plt.title('Indoor Temperature Over Time')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(time_array, T_water_out_list)
plt.xlabel('Time (hours)')
plt.ylabel('Water Outlet Temp (°C)')
plt.title('Chilled Water Outlet Temperature Over Time')
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 3))
plt.plot(time_array, cooling_status_list, drawstyle='steps-post', color='green')
plt.yticks([0, 1], ['OFF', 'ON'])
plt.xlabel('Time (hours)')
plt.ylabel('Cooling System')
plt.title('Cooling System ON/OFF Over Time')
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(time_array, pmv_list, label='PMV', color='orange')
plt.axhline(0.5, color='green', linestyle='--')
plt.axhline(-0.5, color='green', linestyle='--')
plt.xlabel('Time (hours)')
plt.ylabel('PMV')
plt.title('Predicted Mean Vote Over Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 4))
plt.plot(time_array, ppd_list, label='PPD (%)', color='red')
plt.axhline(10, color='blue', linestyle='--', label='Comfort Threshold')
plt.xlabel('Time (hours)')
plt.ylabel('PPD (%)')
plt.title('Predicted Percentage Dissatisfied Over Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
