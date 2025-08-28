# Fick's 2nd Law Nonsteady-State Diffusion Calculator
# Tabulated z and erf(z) values from research paper

erf_table = [
	[0.00,  0.0000],
	[0.025, 0.0282],
	[0.05,  0.0564],
	[0.10,  0.1125],
	[0.15,  0.1680],
	[0.20,  0.2227],
	[0.25,  0.2763],
	[0.30,  0.3286],
	[0.35,  0.3794],
	[0.40,  0.4284],
	[0.45,  0.4755],
	[0.50,  0.5205],
	[0.55,  0.5633],
	[0.60,  0.6039],
	[0.65,  0.6420],
	[0.70,  0.6778],
	[0.75,  0.7112],
	[0.80,  0.7421],
	[0.85,  0.7707],
	[0.90,  0.7970],
	[0.95,  0.8209],
	[1.00,  0.8427],
	[1.10,  0.8802],
	[1.20,  0.9103],
	[1.30,  0.9340],
	[1.40,  0.9523],
	[1.50,  0.9661],
	[1.60,  0.9763],
	[1.70,  0.9838],
	[1.80,  0.9891],
	[1.90,  0.9928],
	[2.00,  0.9953],
	[2.20,  0.9981],
	[2.40,  0.9993],
	[2.60,  0.9998],
	[2.80,  0.9999]
]

import os
import logging
from datetime import datetime
import math

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
log_path = os.path.join(script_dir, 'ficks_2nd_law_log.txt')

# Setup logging to write in the same folder as the script
logging.basicConfig(filename=log_path, level=logging.INFO, format='%(message)s')

def log_block(header, data_dict):
	now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	logging.info("\n" + "="*50)
	logging.info(f"{now} | {header}")
	logging.info("-"*50)
	for key, value in data_dict.items():
		logging.info(f"{key:<30}: {value}")
	logging.info("="*50 + "\n")

def safe_float(prompt):
	while True:
		try:
			return float(input(prompt))
		except ValueError:
			print("Invalid input. Please enter a numeric value.")

def ficks_second_law():
	while True:
		print("\nFick's Second Law Nonsteady-State Diffusion Calculator")
		print("Which variable do you want to solve for? (Cxt, C0, Cs, x, D or t)")
		target = input("Enter variable: ").strip().lower()

		try:
			if target == 'cxt':
				C0 = safe_float("Enter initial concentration C0 [wt%]: ")
				Cs = safe_float("Enter surface concentration Cs [wt%]: ")
				x = safe_float("Enter distance x [m]: ")
				D = safe_float("Enter diffusion coefficient D [m^2/s]: ")
				t = safe_float("Enter time t [s]: ")
				if Cs == C0:
					print("Error: (Cs - C0) cannot be zero (division by zero).")
					continue
				denom = 2 * math.sqrt(D * t)
				if denom == 0:
					print("Error: (2 * sqrt(D * t)) cannot be zero (division by zero).")
					continue
				z = x / denom
				if not (0 <= z <= 2.8):
					print("Error: z must be between 0 and 2.8.")
					continue
				erf_z = round(math.erf(z), 4)
				if not (0 <= erf_z <= 0.9999):
					print("Error: erf(z) must be between 0 and 0.9999.")
					continue
				Cxt = C0 + ((Cs - C0) * (1 - erf_z))
				print(f"z = {z:.4e}")
				print(f"erf(z) = {erf_z:.4f}")
				print(f"Concentration Cxt [wt%] = {Cxt:.4e}")
				log_block("Calculated Concentration Cxt [wt%]", {
					"C0 [wt%]": C0, "Cs [wt%]": Cs, "x [m]": x, "D [m^2/s]": D, "t [s]": t, "z": z, "erf(z)": erf_z, "Cxt [wt%]": Cxt
				})
			elif target == 'c0':
				Cxt = safe_float("Enter concentration Cxt [wt%]: ")
				Cs = safe_float("Enter surface concentration Cs [wt%]: ")
				x = safe_float("Enter distance x [m]: ")
				D = safe_float("Enter diffusion coefficient D [m^2/s]: ")
				t = safe_float("Enter time t [s]: ")
				denom = 2 * math.sqrt(D * t)
				if denom == 0:
					print("Error: (2 * sqrt(D * t)) cannot be zero (division by zero).")
					continue
				z = x / denom
				if not (0 <= z <= 2.8):
					print("Error: z must be between 0 and 2.8.")
					continue
				erf_z = round(math.erf(z), 4)
				if not (0 <= erf_z <= 0.9999):
					print("Error: erf(z) must be between 0 and 0.9999.")
					continue
				a = 1 - erf_z
				C0 = (Cxt - (a * Cs)) / (1 - a) if (1 - a) != 0 else 0
				if Cs == C0:
					print("Error: (Cs - C0) cannot be zero (division by zero).")
					continue
				print(f"z = {z:.4e}")
				print(f"erf(z) = {erf_z:.4f}")
				print(f"Initial concentration C0 [wt%]= {C0:.4e}")
				log_block("Calculated Initial Concentration C0 [wt%]", {
					"Cxt [wt%]": Cxt, "Cs [wt%]": Cs, "x [m]": x, "D [m^2/s]": D, "t [s]": t, "z": z, "erf(z)": erf_z, "C0 [wt%]": C0
				})

			elif target == 'cs':
				Cxt = safe_float("Enter concentration Cxt [wt%]: ")
				C0 = safe_float("Enter initial concentration C0 [wt%]: ")
				x = safe_float("Enter distance x [m]: ")
				D = safe_float("Enter diffusion coefficient D [m^2/s]: ")
				t = safe_float("Enter time t [s]: ")
				denom = 2 * math.sqrt(D * t)
				if denom == 0:
					print("Error: (2 * sqrt(D * t)) cannot be zero (division by zero).")
					continue
				z = x / denom
				if not (0 <= z <= 2.8):
					print("Error: z must be between 0 and 2.8.")
					continue
				erf_z = round(math.erf(z), 4)
				if not (0 <= erf_z <= 0.9999):
					print("Error: erf(z) must be between 0 and 0.9999.")
					continue
				Cs = ((Cxt - C0) / (1 - erf_z)) + C0 if z != 1 else None
				if Cs == C0:
					print("Error: (Cs - C0) cannot be zero (division by zero).")
					continue
				print(f"z = {z:.4e}")
				print(f"erf(z) = {erf_z:.4f}")
				print(f"Surface concentration Cs [wt%]= {Cs:.4e}")
				log_block("Calculated Surface Concentration Cs [wt%]", {
					"Cxt [wt%]": Cxt, "C0 [wt%]": C0, "x [m]": x, "D [m^2/s]": D, "t [s]": t, "z": z, "erf(z)": erf_z, "Cs [wt%]": Cs
				})

			elif target == 'x':
				Cxt = safe_float("Enter concentration Cxt [wt%]: ")
				C0 = safe_float("Enter initial concentration C0 [wt%]: ")
				Cs = safe_float("Enter surface concentration Cs [wt%]: ")
				D = safe_float("Enter diffusion coefficient D [m^2/s]: ")
				t = safe_float("Enter time t [s]: ")
				if Cs == C0:
					print("Error: (Cs - C0) cannot be zero (division by zero).")
					again = input("\nWould you like to use the calculator again? (y/n): ").strip().lower()
					if again != 'y':
						print("Thank you for using the Fick's Second Law calculator!")
						return
					else:
						continue
				cerfz = 1 - ((Cxt - C0) / (Cs - C0))
				# Find closest lower and higher values in erf_table
				lower = None
				higher = None
				for i in range(len(erf_table)):
					if erf_table[i][1] <= cerfz:
						lower = erf_table[i]
					if erf_table[i][1] >= cerfz and higher is None:
						higher = erf_table[i]
				if lower is None or higher is None:
					print("Calculated erf(z) is out of table range.")
					again = input("\nWould you like to use the calculator again? (y/n): ").strip().lower()
					if again != 'y':
						print("Thank you for using the Fick's Second Law calculator!")
						return
					continue
				lz, lerfz = lower
				hz, herfz = higher
				nz = (((cerfz - lerfz) / (herfz - lerfz)) * (hz - lz)) + lz if herfz != lerfz else lz
				x = nz * 2 * math.sqrt(D * t)
				print(f"Calculated z = {cerfz:.4f}")
				print(f"Lower: z = {lz}, erf(z) = {lerfz}")
				print(f"Higher: z = {hz}, erf(z) = {herfz}")
				print(f"Interpolated z = {nz:.4f}")
				print(f"D [m^2/s] = {D:.4e}")
				print(f"t [s] = {t:.4e}")
				print(f"Calculated x [m] = {x:.4e}")
				log_block("Calculated Distance x [m]", {
					"erf(z)": cerfz, "z": nz, "D [m^2/s]": D, "t [s]": t, "x [m]": x
				})
			elif target == 't':
				Cxt = safe_float("Enter concentration Cxt [wt%]: ")
				C0 = safe_float("Enter initial concentration C0 [wt%]: ")
				Cs = safe_float("Enter surface concentration Cs [wt%]: ")
				D = safe_float("Enter diffusion coefficient D [m^2/s]: ")
				x = safe_float("Enter distance x [m]: ")
				if Cs == C0:
					print("Error: (Cs - C0) cannot be zero (division by zero).")
					return
				cerfz = 1 - ((Cxt - C0) / (Cs - C0))
				lower = None
				higher = None
				for i in range(len(erf_table)):
					if erf_table[i][1] <= cerfz:
						lower = erf_table[i]
					if erf_table[i][1] >= cerfz and higher is None:
						higher = erf_table[i]
				if lower is None or higher is None:
					print("Calculated erf(z) is out of table range.")
					return
				lz, lerfz = lower
				hz, herfz = higher
				nz = (((cerfz - lerfz) / (herfz - lerfz)) * (hz - lz)) + lz if herfz != lerfz else lz
				t = (x / (2 * nz * math.sqrt(D))) ** 2 if nz != 0 and D > 0 else 0
				if t < 0:
					print("Calculated t is negative, check your inputs.")
					again = input("\nWould you like to use the calculator again? (y/n): ").strip().lower()
					if again != 'y':
						print("Thank you for using the Fick's Second Law calculator!")
						break
					continue
				print(f"Calculated erf(z) = {cerfz:.4f}")
				print(f"Lower: z = {lz}, erf(z) = {lerfz}")
				print(f"Higher: z = {hz}, erf(z) = {herfz}")
				print(f"Interpolated z = {nz:.4f}")
				print(f"D [m^2/s] = {D:.4e}")
				print(f"x [m] = {x:.4e}")
				print(f"Calculated t [s] = {t:.4e}")
				log_block("Calculated Time t [s]", {
					"erf(z)": cerfz, "z": nz, "D [m^2/s]": D, "x [m]": x, "t [s]": t
				})
			elif target == 'd':
				Cxt = safe_float("Enter concentration Cxt [wt%]: ")
				C0 = safe_float("Enter initial concentration C0 [wt%]: ")
				Cs = safe_float("Enter surface concentration Cs [wt%]: ")
				t = safe_float("Enter time t [s]: ")
				x = safe_float("Enter distance x [m]: ")
				if Cs == C0:
					print("Error: (Cs - C0) cannot be zero (division by zero).")
					return
				cerfz = 1 - ((Cxt - C0) / (Cs - C0))
				lower = None
				higher = None
				for i in range(len(erf_table)):
					if erf_table[i][1] <= cerfz:
						lower = erf_table[i]
					if erf_table[i][1] >= cerfz and higher is None:
						higher = erf_table[i]
				if lower is None or higher is None:
					print("Calculated erf(z) is out of table range.")
					return
				lz, lerfz = lower
				hz, herfz = higher
				nz = (((cerfz - lerfz) / (herfz - lerfz)) * (hz - lz)) + lz if herfz != lerfz else lz
				D = (x / (2 * nz * math.sqrt(t))) ** 2 if nz != 0 and t > 0 else 0
				if D < 0:
					print("Calculated D is negative, check your inputs.")
					again = input("\nWould you like to use the calculator again? (y/n): ").strip().lower()
					if again != 'y':
						print("Thank you for using the Fick's Second Law calculator!")
						break
					continue
				print(f"Calculated erf(z) = {cerfz:.4f}")
				print(f"Lower: z = {lz}, erf(z) = {lerfz}")
				print(f"Higher: z = {hz}, erf(z) = {herfz}")
				print(f"Interpolated z = {nz:.4f}")
				print(f"t [s] = {t:.4e}")
				print(f"x [m] = {x:.4e}")
				print(f"Calculated D [m^2/s] = {D:.4e}")
				log_block("Calculated Diffusion Coefficient D [m^2/s]", {
					"erf(z)": cerfz, "z": nz, "t [s]": t, "x [m]": x, "D [m^2/s]": D
				})
			else:
				print("Invalid selection. Please choose Cxt, C0, Cs, x, D or t.")

		except Exception as e:
			print(f"An unexpected error occurred: {e}")
			again = input("\nWould you like to use the calculator again? (y/n): ").strip().lower()
			if again != 'y':
				print("Thank you for using the Fick's Second Law calculator!")
				break
			else:
				continue

		again = input("\nWould you like to use the calculator again? (y/n): ").strip().lower()
		if again != 'y':
			print("Thank you for using the Fick's Second Law calculator!")
			break

if __name__ == "__main__":
	ficks_second_law()
