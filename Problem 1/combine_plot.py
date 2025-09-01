import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import os

# --- Configuration ---
KINETICS_SUMMARY_CSV = 'final_integral_kinetics_summary.csv'
GAS_CONSTANT_R = 8.314  # J/(mol·K)

def generate_combined_arrhenius_plots():
    """
    Reads the final kinetics summary and for each reaction (RX1, RX2, RX3),
    generates a single, combined Arrhenius plot showing the data and model fits
    (linear and polynomial) for both k1 and k2.
    """
    print("🚀 Starting Combined & Comparative Arrhenius Analysis...")

    try:
        df_kinetics = pd.read_csv(KINETICS_SUMMARY_CSV)
    except FileNotFoundError:
        print(f"❌ ERROR: '{KINETICS_SUMMARY_CSV}' not found.")
        return

    # --- Main Loop: Iterate through each Reaction ID ---
    for reaction_id in df_kinetics['Reaction_ID'].unique():
        output_folder = f'{reaction_id}_final_arrhenius_analysis'
        os.makedirs(output_folder, exist_ok=True)
        print(f"\n{'='*50}\nProcessing Reaction: {reaction_id}\n{'='*50}")

        df_rx = df_kinetics[df_kinetics['Reaction_ID'] == reaction_id]

        # --- 1. Prepare Data for k1 and k2 ---
        df_agg = df_rx.groupby('Temperature_C')[['k1', 'k2']].mean().reset_index()
        df_agg = df_agg[(df_agg['k1'] > 0) & (df_agg['k2'] > 0)] # Ensure k values are positive
        
        if len(df_agg) < 3:
            print("  - Skipping: Not enough temperature data points for a meaningful fit.")
            continue

        df_agg['Temp_K'] = df_agg['Temperature_C'] + 273.15
        df_agg['1/T'] = 1 / df_agg['Temp_K']
        df_agg['ln(k1)'] = np.log(df_agg['k1'])
        df_agg['ln(k2)'] = np.log(df_agg['k2'])

        # --- 2. Fit Models for k1 ---
        x_data, y_data_k1 = df_agg['1/T'], df_agg['ln(k1)']
        slope_lin1, intercept_lin1, r_val_lin1, _, _ = linregress(x_data, y_data_k1)
        Ea1_lin = -slope_lin1 * GAS_CONSTANT_R / 1000
        poly_coeffs1 = np.polyfit(x_data, y_data_k1, 2)
        poly_fit_func1 = np.poly1d(poly_coeffs1)
        r_sq_poly1 = 1 - (np.sum((y_data_k1 - poly_fit_func1(x_data))**2) / np.sum((y_data_k1 - np.mean(y_data_k1))**2))
        
        # --- 3. Fit Models for k2 ---
        y_data_k2 = df_agg['ln(k2)']
        slope_lin2, intercept_lin2, r_val_lin2, _, _ = linregress(x_data, y_data_k2)
        Ea2_lin = -slope_lin2 * GAS_CONSTANT_R / 1000
        poly_coeffs2 = np.polyfit(x_data, y_data_k2, 2)
        poly_fit_func2 = np.poly1d(poly_coeffs2)
        r_sq_poly2 = 1 - (np.sum((y_data_k2 - poly_fit_func2(x_data))**2) / np.sum((y_data_k2 - np.mean(y_data_k2))**2))

        print(f"  - k1 Fits: Linear R²={r_val_lin1**2:.4f}, Polynomial R²={r_sq_poly1:.4f}")
        print(f"  - k2 Fits: Linear R²={r_val_lin2**2:.4f}, Polynomial R²={r_sq_poly2:.4f}")

        # --- 4. Generate Combined Plot ---
        plt.figure(figsize=(14, 9))
        x_fine = np.linspace(x_data.min(), x_data.max(), 200)

        # Plot k1 data and fits
        plt.scatter(x_data, y_data_k1, color='blue', marker='o', s=100, label='ln(k₁) Data')
        plt.plot(x_fine, slope_lin1 * x_fine + intercept_lin1, 'b-', linewidth=2, label='k₁ Standard Model (Linear)')
        plt.plot(x_fine, poly_fit_func1(x_fine), 'b--', linewidth=2, label='k₁ Temp-Dependent Model (Poly)')

        # Plot k2 data and fits
        plt.scatter(x_data, y_data_k2, color='red', marker='s', s=100, label='ln(k₂) Data')
        plt.plot(x_fine, slope_lin2 * x_fine + intercept_lin2, 'r-', linewidth=2, label='k₂ Standard Model (Linear)')
        plt.plot(x_fine, poly_fit_func2(x_fine), 'r--', linewidth=2, label='k₂ Temp-Dependent Model (Poly)')
        
        plt.title(f'Combined Arrhenius Analysis for {reaction_id}', fontsize=18)
        plt.xlabel('1 / Temperature (K⁻¹)', fontsize=14)
        plt.ylabel('ln(k)', fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(True, linestyle='--')

        # --- Add a detailed summary text box ---
        info_text = (f'--- Standard Model (Linear Fit) ---\n'
                     f'$E_{{a1}}$ = {Ea1_lin:.2f} kJ/mol (R² = {r_val_lin1**2:.3f})\n'
                     f'$E_{{a2}}$ = {Ea2_lin:.2f} kJ/mol (R² = {r_val_lin2**2:.3f})\n\n'
                     f'--- Temp-Dependent Model (Poly Fit) ---\n'
                     f'k₁ Fit R² = {r_sq_poly1:.3f}\n'
                     f'k₂ Fit R² = {r_sq_poly2:.3f}')
        plt.text(0.03, 0.97, info_text, transform=plt.gca().transAxes, fontsize=11,
                 va='top', bbox=dict(boxstyle='round', fc='wheat', alpha=0.8))

        filename = f'{reaction_id}_combined_arrhenius_plot.png'
        filepath = os.path.join(output_folder, filename)
        plt.savefig(filepath)
        plt.close()

    print(f"\n\n{'='*50}\n✅ COMBINED ARRHENIUS ANALYSIS COMPLETE\n{'='*50}")
    print("📈 All combined plots have been saved into their respective folders.")

if __name__ == '__main__':
    generate_combined_arrhenius_plots()
