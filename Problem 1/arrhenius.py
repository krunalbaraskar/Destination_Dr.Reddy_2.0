import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os

# --- Configuration ---
KINETICS_SUMMARY_CSV = 'final_integral_kinetics_summary.csv'
GAS_CONSTANT_R = 8.314  # J/(mol·K)

def generate_advanced_arrhenius_plots():
    """
    Reads the final kinetics summary, and for each reaction (RX1, RX2, RX3) and
    each rate constant (k1, k2), it generates an Arrhenius plot comparing a
    standard linear fit with a temperature-dependent (2nd order polynomial) fit.
    """
    print("🚀 Starting Advanced Arrhenius Analysis for all reactions...")

    try:
        df_kinetics = pd.read_csv(KINETICS_SUMMARY_CSV)
    except FileNotFoundError:
        print(f"❌ ERROR: '{KINETICS_SUMMARY_CSV}' not found.")
        return

    # --- Main Loop: Iterate through each Reaction ID ---
    for reaction_id in df_kinetics['Reaction_ID'].unique():
        output_folder = f'{reaction_id}_advanced_arrhenius_plots'
        os.makedirs(output_folder, exist_ok=True)
        print(f"\n{'='*50}\nProcessing Reaction: {reaction_id}\n{'='*50}")

        df_rx = df_kinetics[df_kinetics['Reaction_ID'] == reaction_id]

        # Process both k1 and k2
        for k_num in [1, 2]:
            k_col = f'k{k_num}'
            
            # --- 1. Prepare Data ---
            # Aggregate data to get one point per temperature
            df_agg = df_rx.groupby('Temperature_C')[k_col].mean().reset_index()
            
            # Filter out any non-positive k values before taking log
            df_agg = df_agg[df_agg[k_col] > 0]
            if len(df_agg) < 3:
                print(f"  - Skipping k{k_num}: Not enough data points for a meaningful fit.")
                continue

            df_agg['Temp_K'] = df_agg['Temperature_C'] + 273.15
            df_agg['1/T'] = 1 / df_agg['Temp_K']
            df_agg[f'ln(k{k_num})'] = np.log(df_agg[k_col])

            x_data = df_agg['1/T']
            y_data = df_agg[f'ln(k{k_num})']
            
            # --- 2. Fit Both Models ---
            # Model 1: Standard Arrhenius (Linear Fit)
            slope_lin, intercept_lin, r_val_lin, _, _ = linregress(x_data, y_data)
            Ea_lin = -slope_lin * GAS_CONSTANT_R / 1000  # in kJ/mol
            y_fit_lin = slope_lin * x_data + intercept_lin

            # Model 2: Temperature-Dependent (2nd Order Polynomial Fit)
            poly_coeffs = np.polyfit(x_data, y_data, 2)
            poly_fit_func = np.poly1d(poly_coeffs)
            y_fit_poly = poly_fit_func(x_data)
            # Calculate R-squared for the polynomial fit
            ss_res_poly = np.sum((y_data - y_fit_poly)**2)
            ss_tot = np.sum((y_data - np.mean(y_data))**2)
            r_sq_poly = 1 - (ss_res_poly / ss_tot)
            
            print(f"  - Analyzing k{k_num}: Linear R²={r_val_lin**2:.4f}, Polynomial R²={r_sq_poly:.4f}")

            # --- 3. Generate Plot ---
            plt.figure(figsize=(12, 8))
            plt.scatter(x_data, y_data, color='black', s=100, zorder=5, label='Experimental Data')
            
            # Plot the fits using a finer x-axis for smooth curves
            x_fine = np.linspace(x_data.min(), x_data.max(), 200)
            plt.plot(x_fine, slope_lin * x_fine + intercept_lin, 'b-', linewidth=2,
                     label=f'Standard Model (Linear Fit, R²={r_val_lin**2:.3f})')
            plt.plot(x_fine, poly_fit_func(x_fine), 'r--', linewidth=2,
                     label=f'Temp-Dependent Model (Poly Fit, R²={r_sq_poly:.3f})')
            
            plt.title(f'Arrhenius Model Comparison for {reaction_id} - k{k_num}', fontsize=16)
            plt.xlabel('1 / Temperature (K⁻¹)', fontsize=12)
            plt.ylabel(f'ln(k{k_num})', fontsize=12)
            plt.legend()
            plt.grid(True, linestyle='--')

            info_text = (f'Standard Model (Linear):\n'
                         f'$E_a$ = {Ea_lin:.2f} kJ/mol')
            plt.text(0.05, 0.95, info_text, transform=plt.gca().transAxes, fontsize=12,
                     va='top', bbox=dict(boxstyle='round', fc='wheat', alpha=0.8))

            filename = f'{reaction_id}_arrhenius_comparison_k{k_num}.png'
            filepath = os.path.join(output_folder, filename)
            plt.savefig(filepath)
            plt.close()

    print(f"\n\n{'='*50}\n✅ ADVANCED ARRHENIUS ANALYSIS COMPLETE\n{'='*50}")
    print("📈 All comparison plots have been saved into their respective folders.")

if __name__ == '__main__':
    generate_advanced_arrhenius_plots()
