import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import os

# --- Configuration ---
REACTION_FILES = ['RX1.csv', 'RX2.csv', 'RX3.csv']
KINETICS_SUMMARY_CSV = 'final_integral_kinetics_summary.csv'
FRESH_FEED_FLOW_RATE_Q_F = 1000 # ml/h for production rate calculation

def generate_final_graphs():
    """
    Reads a summary of kinetic parameters (n1, n2, k1, k2, type) and generates
    final, polished plots for each experiment, showing the model fit and all
    key performance metrics.
    """
    print("🚀 Starting Final Report Graph Generation...")

    try:
        df_kinetics = pd.read_csv(KINETICS_SUMMARY_CSV)
        df_raw_map = {os.path.splitext(f)[0]: pd.read_csv(f) for f in REACTION_FILES}
    except FileNotFoundError as e:
        print(f"❌ ERROR: Missing a required file: {e}. Please ensure all necessary CSV files are present.")
        return

    # --- Calculate Average Orders for each Reaction ID ---
    avg_orders = df_kinetics.groupby('Reaction_ID')[['n1', 'n2']].mean().to_dict('index')

    # --- Main Loop ---
    for _, row in df_kinetics.iterrows():
        # 1. Extract parameters for the current experiment
        reaction_id, temp, conc, r_type = row['Reaction_ID'], row['Temperature_C'], row['Initial_Conc_A'], row['Model_Type']
        k1, k2 = row['k1'], row['k2']
        
        # Use the average orders for the model
        avg_n1 = avg_orders[reaction_id]['n1']
        avg_n2 = avg_orders[reaction_id]['n2']
        
        print(f"\nPlotting: {reaction_id} @ {temp}°C, C_A0={conc} (Model: {r_type}, n1={avg_n1:.2f}, n2={avg_n2:.2f})")
        
        output_folder = f'{reaction_id}_final_report_graphs'
        os.makedirs(output_folder, exist_ok=True)
        
        subset_df = df_raw_map[reaction_id][(df_raw_map[reaction_id]['Temperature'] == temp) & (df_raw_map[reaction_id]['Initial Concentration'] == conc)]

        # --- 2. Define the correct ODE model based on type ---
        if r_type == 'Series':
            def model_ode(C, t, k1, n1, k2, n2):
                C_A, C_B, _ = C
                dC_A_dt = -k1 * (max(C_A, 1e-9)**n1)
                dC_B_dt = (k1 * (max(C_A, 1e-9)**n1)) - (k2 * (max(C_B, 1e-9)**n2))
                dC_I_dt = k2 * (max(C_B, 1e-9)**n2)
                return [dC_A_dt, dC_B_dt, dC_I_dt]
        else: # Parallel
            def model_ode(C, t, k1, n1, k2, n2):
                C_A, _, _ = C
                rate_B = k1 * (max(C_A, 1e-9)**n1)
                rate_I = k2 * (max(C_A, 1e-9)**n2)
                dC_A_dt = -rate_B - rate_I; dC_B_dt = rate_B; dC_I_dt = rate_I
                return [dC_A_dt, dC_B_dt, dC_I_dt]

        # --- 3. Simulate and Analyze the Model Curve ---
        time_fine = np.linspace(0, subset_df['Time'].max(), 500)
        model_curve = odeint(model_ode, [conc, 0, 0], time_fine, args=(k1, avg_n1, k2, avg_n2))
        
        # Find B_max from the smooth model curve
        b_max_idx = np.argmax(model_curve[:, 1])
        time_at_b_max = time_fine[b_max_idx]
        conc_b_max = model_curve[b_max_idx, 1]
        conc_i_at_b_max = model_curve[b_max_idx, 2]

        # Calculate metrics from the model's peak
        yield_at_b_max = (conc_b_max / conc) * 100.0
        selectivity_at_b_max = conc_b_max / conc_i_at_b_max if conc_i_at_b_max > 1e-6 else float('inf')
        
        # Calculate overall metrics from experimental data
        final_row = subset_df.iloc[-1]
        overall_yield = (final_row['B'] / conc) * 100.0
        overall_prod_rate = final_row['B'] * FRESH_FEED_FLOW_RATE_Q_F

        # --- 4. Generate Final Plot ---
        plt.figure(figsize=(12, 8))
        plt.plot(subset_df['Time'], subset_df['A'], 'o', color='blue', label='Reactant [A]')
        plt.plot(subset_df['Time'], subset_df['B'], 's', color='green', label='Product [B]')
        plt.plot(subset_df['Time'], subset_df['I'], '^', color='red', label='Impurity [I]')
        
        plt.plot(time_fine, model_curve[:, 0], 'b-')
        plt.plot(time_fine, model_curve[:, 1], 'g-')
        plt.plot(time_fine, model_curve[:, 2], 'r-')

        plt.axvline(x=time_at_b_max, color='gray', linestyle='--', label=f'Time for Max [B] = {time_at_b_max:.2f} h')
        plt.plot(time_at_b_max, conc_b_max, '*', markersize=18, color='gold', markeredgecolor='black', label=f'Max [B] = {conc_b_max:.2f}')

        # Add overall metrics to the legend
        plt.plot([], [], ' ', label='--------------------')
        plt.plot([], [], ' ', label=f'Overall Yield: {overall_yield:.1f}%')
        plt.plot([], [], ' ', label=f'Overall B Produced: {overall_prod_rate:.0f} mg/h')
        
        plt.title(f'{r_type} Reaction: {reaction_id} at {temp}°C (Initial [A] = {conc} mg/ml)', fontsize=16)
        plt.xlabel('Time (h)'); plt.ylabel('Concentration (mg/ml)')
        plt.legend(loc='best'); plt.grid(True, linestyle='--')

        info_text = (f'Model: {r_type} (n₁={avg_n1:.2f}, n₂={avg_n2:.2f})\n\n'
                     f'$k_1$ = {k1:.4f} \n'
                     f'$k_2$ = {k2:.4f} \n\n'
                     f'At Model $B_{{max}}$ (t={time_at_b_max:.2f}h):\n'
                     f'Yield = {yield_at_b_max:.1f}%\n'
                     f'Selectivity = {selectivity_at_b_max:.1f}')
        plt.text(0.95, 0.65, info_text, transform=plt.gca().transAxes, ha='right', va='top',
                 bbox=dict(boxstyle='round', fc='wheat', alpha=0.8))

        filename = f'{reaction_id}_T{temp}_C{conc}_final_model.png'
        filepath = os.path.join(output_folder, filename)
        plt.savefig(filepath)
        plt.close()

    print(f"\n\n{'='*50}\n✅ FINAL REPORT GRAPH GENERATION COMPLETE\n{'='*50}")
    print("📈 All final model-based plots have been saved to their respective folders.")

if __name__ == '__main__':
    generate_final_graphs()
