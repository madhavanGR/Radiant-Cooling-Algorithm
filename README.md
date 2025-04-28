
# Thermal Simulation of Wall Cooling with Chilled Water Pipes

This project simulates the thermal behavior of a vertical wall embedded with chilled water pipes to control the indoor temperature of a room. It also calculates human thermal comfort metrics (PMV and PPD) over a 24-hour period based on external weather conditions.

## Project Structure

- `thermal_simulation.py`: Main Python script that runs the simulation.
- `bengaluru_hourly_temps.csv`: External hourly weather data (Temperature and Humidity for Bengaluru).

## Key Features

- Finite difference simulation of wall heat conduction
- Chilled water loop for cooling
- Indoor air thermal dynamics
- HVAC control (On/Off based on indoor setpoints)
- Comfort assessment using PMV and PPD indices
- Visualization of temperature profiles, cooling system status, and comfort metrics

## Simulation Details

- **Wall size**: 0.2 m thickness × 2.5 m height
- **Wall material**: Thermal conductivity = 0.8 W/mK
- **Pipe properties**: U-value = 200 W/m²K, pipe spacing = 0.2 m
- **Indoor air control**:
  - Setpoint High: 24°C
  - Setpoint Low: 22°C
- **Time step**: 60 seconds
- **Simulation duration**: 24 hours

## How to Run

1. Install dependencies:
```bash
pip install numpy pandas matplotlib
```

2. Make sure `bengaluru_hourly_temps.csv` is placed in the same directory as the script.

3. Run the simulation:
```bash
python thermal_simulation.py
```

4. The script will generate multiple plots:
- Indoor temperature vs. time
- Chilled water outlet temperature vs. time
- Cooling system ON/OFF status
- PMV (Predicted Mean Vote)
- PPD (Predicted Percentage Dissatisfied)

## Inputs Required

- **External Weather File**: `bengaluru_hourly_temps.csv`
  - Should contain two columns: `Temperature` and `Humidity`
  - 24 rows representing each hour

## Future Improvements

- Add multi-wall/multi-zone modeling
- Integrate energy consumption tracking (kWh)
- Simulate variable flow rates and dynamic cooling
- Include real-time animation of wall temperature fields

## License

This project is released under the MIT License.

---

✨ Developed to demonstrate building thermal performance and comfort modeling using Python!
