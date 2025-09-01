import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
import os

# --- Configuration ---
REACTION_FILES = ['RX1.csv', 'RX2.csv', 'RX3.csv']
TYPE_SUMMARY_CSV = 'reaction_type_summary.csv'
FINAL_KINETICS_SUMMARY_CSV = 'final_integral_kinetics_summary.csv'

# --- Integral Analysis Functions ---
def integral_series(subset_df):
    """Performs global optimization for a Series model (A->B->I)."""
    time = subset_df['Time'].values
    exp_data = subset_df[['A', 'B', 'I']].values
    C_A0 = subset_df['Initial Concentration'].iloc[0]
    
    def model(C, t, k1, n1, k2, n2):
        C_A, C_B, _ = C
        dC_A_dt = -k1 * (max(C_A, 1e-9)**n1)
        dC_B_dt = (k1 * (max(C_A, 1e-9)**n1)) - (k2 * (max(C_B, 1e-9)**n2))
        dC_I_dt = k2 * (max(C_B, 1e-9)**n2)
        return [dC_A_dt, dC_B_dt, dC_I_dt]

    def objective(params):
        k1, n1, k2, n2 = params
        pred = odeint(model, [C_A0, 0, 0], time, args=(k1, n1, k2, n2))
        return np.sum((exp_data - pred)**2)

    # Initial guesses [k1, n1, k2, n2] and bounds
    initial_guesses = [1.0, 1.0, 0.1, 1.0]
    bounds = [(1e-6, None), (0, 3), (1e-6, None), (0, 3)]
    res = minimize(objective, initial_guesses, method='L-BFGS-B', bounds=bounds)
    
    k1_fit, n1_fit, k2_fit, n2_fit = res.x
    return {'n1': n1_fit, 'k1': k1_fit, 'n2': n2_fit, 'k2': k2_fit, 'ssr': res.fun, 'model_func': model}

def integral_parallel(subset_df):
    """Performs global optimization for a Parallel model (A->B, A->I)."""
    time = subset_df['Time'].values
    exp_data = subset_df[['A', 'B', 'I']].values
    C_A0 = subset_df['Initial Concentration'].iloc[0]

    def model(C, t, k1, n1, k2, n2):
        C_A, _, _ = C
        rate_B = k1 * (max(C_A, 1e-9)**n1)
        rate_I = k2 * (max(C_A, 1e-9)**n2)
        dC_A_dt = -rate_B - rate_I
        dC_B_dt = rate_B
        dC_I_dt = rate_I
        return [dC_A_dt, dC_B_dt, dC_I_dt]

    def objective(params):
        k1, n1, k2, n2 = params
        pred = odeint(model, [C_A0, 0, 0], time, args=(k1, n1, k2, n2))
        return np.sum((exp_data - pred)**2)

    initial_guesses = [1.0, 1.0, 0.1, 1.0]
    bounds = [(1e-6, None), (0, 3), (1e-6, None), (0, 3)]
    res = minimize(objective, initial_guesses, method='L-BFGS-B', bounds=bounds)
    
    k1_fit, n1_fit, k2_fit, n2_fit = res.x
    return {'n1': n1_fit, 'k1': k1_fit, 'n2': n2_fit, 'k2': k2_fit, 'ssr': res.fun, 'model_func': model}


def main_analyzer():
    print("🚀 Starting Comprehensive Integral Kinetic Analysis...")
    try:
        df_type = pd.read_csv(TYPE_SUMMARY_CSV)
        df_raw_map = {os.path.splitext(f)[0]: pd.read_csv(f) for f in REACTION_FILES}
    except FileNotFoundError as e:
        print(f"❌ ERROR: Missing a required file: {e}. Please ensure all CSV files are present.")
        return

    final_results = []
    
    for _, row in df_type.iterrows():
        reaction_id, temp, conc, r_type = row['Reaction_ID'], row['Temperature_C'], row['Initial_Conc_A_mg/ml'], row['Determined_Reaction_Type']
        print(f"\nAnalyzing: {reaction_id} @ {temp}°C, C_A0={conc} (Model: {r_type})")
        
        output_folder = f'{reaction_id}_final_integral_fits'
        os.makedirs(output_folder, exist_ok=True)
        
        subset_df = df_raw_map[reaction_id][(df_raw_map[reaction_id]['Temperature'] == temp) & (df_raw_map[reaction_id]['Initial Concentration'] == conc)]
        
        if r_type == 'Series':
            int_res = integral_series(subset_df)
        else: # Parallel
            int_res = integral_parallel(subset_df)

        # Combine results and save
        k1, n1, k2, n2 = int_res['k1'], int_res['n1'], int_res['k2'], int_res['n2']
        full_res = {
            'Reaction_ID': reaction_id, 'Temperature_C': temp, 'Initial_Conc_A': conc,
            'Model_Type': r_type, 'n1': n1, 'k1': k1, 'n2': n2, 'k2': k2, 'SSR': int_res['ssr']
        }
        final_results.append(full_res)

        # Plot Integral Fit
        time_fine = np.linspace(0, subset_df['Time'].max(), 200)
        pred_curve = odeint(int_res['model_func'], [conc, 0, 0], time_fine, args=(k1, n1, k2, n2))

        plt.figure(figsize=(12, 8))
        plt.scatter(subset_df['Time'], subset_df['A'], c='blue', label='Exp. [A]')
        plt.scatter(subset_df['Time'], subset_df['B'], c='green', label='Exp. [B]')
        plt.scatter(subset_df['Time'], subset_df['I'], c='red', label='Exp. [I]')
        plt.plot(time_fine, pred_curve[:,0], 'b-', label='Model Fit [A]')
        plt.plot(time_fine, pred_curve[:,1], 'g-', label='Model Fit [B]')
        plt.plot(time_fine, pred_curve[:,2], 'r-', label='Model Fit [I]')
        
        plt.title(f'Integral Model Fit: {reaction_id} @ {temp}°C (C_A0={conc}) - {r_type} Model', fontsize=16)
        plt.xlabel('Time (h)'); plt.ylabel('Concentration (mg/ml)')
        info_text = (f'--- Final Fit Parameters ---\n'
                     f'Step 1 (n₁): {n1:.2f}, k₁: {k1:.4f}\n'
                     f'Step 2 (n₂): {n2:.2f}, k₂: {k2:.4f}\n'
                     f'SSR: {int_res["ssr"]:.2f}')
        plt.text(0.95, 0.95, info_text, transform=plt.gca().transAxes, ha='right', va='top', bbox=dict(fc='wheat', alpha=0.8))
        plt.legend(); plt.grid(True, linestyle='--')
        
        filename = f'{reaction_id}_T{temp}_C{conc}_integral_fit.png'
        filepath = os.path.join(output_folder, filename)
        plt.savefig(filepath)
        plt.close()

    pd.DataFrame(final_results).to_csv(FINAL_KINETICS_SUMMARY_CSV, index=False, float_format='%.4f')
    print(f"\n\n{'='*50}\n✅ FULL KINETIC ANALYSIS COMPLETE\n{'='*50}")
    print(f"📊 A final summary of all kinetic parameters saved to '{FINAL_KINETICS_SUMMARY_CSV}'")

if __name__ == '__main__':
    main_analyzer()
