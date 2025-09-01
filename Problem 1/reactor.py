import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve
import os

# --- Configuration ---
KINETICS_SUMMARY_CSV = 'final_integral_kinetics_summary.csv'
FINAL_COMPARISON_SUMMARY = 'reactor_choice_summary.csv'

# --- PFR Model (ODE Solver) ---
def pfr_model(C, t, k1, n1, k2, n2, model_type):
    C_A, C_B, _ = C
    C_A_safe = max(C_A, 1e-9); C_B_safe = max(C_B, 1e-9)
    if model_type == 'Series':
        dC_A_dt = -k1 * (C_A_safe**n1)
        dC_B_dt = (k1 * (C_A_safe**n1)) - (k2 * (C_B_safe**n2))
        dC_I_dt = k2 * (C_B_safe**n2)
    else: # Parallel
        rate_B = k1 * (C_A_safe**n1); rate_I = k2 * (C_A_safe**n2)
        dC_A_dt = -rate_B - rate_I; dC_B_dt = rate_B; dC_I_dt = rate_I
    return [dC_A_dt, dC_B_dt, dC_I_dt]

# --- CSTR Model (Algebraic Solver) ---
def cstr_model(C, tau, C_A0, k1, n1, k2, n2, model_type):
    C_A, C_B, C_I = C
    if model_type == 'Series':
        eq1 = C_A0 - C_A - (k1 * (C_A**n1) * tau)
        eq2 = 0 - C_B + (k1 * (C_A**n1) * tau) - (k2 * (C_B**n2) * tau)
        eq3 = 0 - C_I + (k2 * (C_B**n2) * tau)
    else: # Parallel
        rate_B = k1 * (C_A**n1); rate_I = k2 * (C_A**n2)
        eq1 = C_A0 - C_A - (rate_B * tau) - (rate_I * tau)
        eq2 = 0 - C_B + (rate_B * tau)
        eq3 = 0 - C_I + (rate_I * tau)
    return [eq1, eq2, eq3]

def compare_reactors():
    print("🚀 Starting PFR vs. CSTR Performance Comparison...")
    try:
        df_kinetics = pd.read_csv(KINETICS_SUMMARY_CSV)
    except FileNotFoundError:
        print(f"❌ ERROR: '{KINETICS_SUMMARY_CSV}' not found.")
        return

    avg_orders = df_kinetics.groupby('Reaction_ID')[['n1', 'n2']].mean().to_dict('index')
    comparison_results = []

    for _, row in df_kinetics.iterrows():
        reaction_id, temp, conc, r_type, k1, k2 = row['Reaction_ID'], row['Temperature_C'], row['Initial_Conc_A'], row['Model_Type'], row['k1'], row['k2']
        n1, n2 = avg_orders[reaction_id]['n1'], avg_orders[reaction_id]['n2']
        
        print(f"\nComparing Reactors for: {reaction_id} @ {temp}°C, C_A0={conc}")

        output_folder = f'{reaction_id}_reactor_comparison_plots'
        os.makedirs(output_folder, exist_ok=True)
        
        res_times = np.linspace(1e-6, 25, 200)

        # --- PFR Simulation ---
        pfr_profile = odeint(pfr_model, [conc, 0, 0], res_times, args=(k1, n1, k2, n2, r_type))
        pfr_yield = (pfr_profile[:, 1] / conc) * 100.0
        pfr_selectivity = pfr_profile[:, 1] / (pfr_profile[:, 2] + 1e-9)

        # --- CSTR Simulation ---
        cstr_profile = np.array([fsolve(cstr_model, x0=[conc*0.5, conc*0.1, conc*0.1], args=(tau, conc, k1, n1, k2, n2, r_type)) for tau in res_times])
        cstr_yield = (cstr_profile[:, 1] / conc) * 100.0
        cstr_selectivity = cstr_profile[:, 1] / (cstr_profile[:, 2] + 1e-9)

        # --- Analyze Results ---
        pfr_max_yield = np.max(pfr_yield)
        cstr_max_yield = np.max(cstr_yield)
        recommended_reactor = 'PFR' if pfr_max_yield > cstr_max_yield else 'CSTR'

        comparison_results.append({
            'Reaction_ID': reaction_id, 'Temperature_C': temp, 'Initial_Conc_A': conc,
            'Model_Type': r_type, 'PFR_Max_Yield_%': pfr_max_yield, 'CSTR_Max_Yield_%': cstr_max_yield,
            'Recommended_Reactor': recommended_reactor
        })
        
        # --- Generate Plots ---
        # Yield Plot
        plt.figure(figsize=(10, 7))
        plt.plot(res_times, pfr_yield, 'b-', label=f'PFR (Max Yield: {pfr_max_yield:.1f}%)')
        plt.plot(res_times, cstr_yield, 'g--', label=f'CSTR (Max Yield: {cstr_max_yield:.1f}%)')
        plt.title(f'Yield Comparison: {reaction_id} @ {temp}°C, C_A0={conc}')
        plt.xlabel('Residence Time (τ) (hours)'); plt.ylabel('Yield of B (%)')
        plt.legend(); plt.grid(True, linestyle='--'); plt.ylim(bottom=0)
        plt.savefig(os.path.join(output_folder, f'{reaction_id}_T{temp}_C{conc}_yield_comparison.png'))
        plt.close()

        # Selectivity Plot
        plt.figure(figsize=(10, 7))
        plt.plot(res_times, pfr_selectivity, 'b-', label='PFR Selectivity')
        plt.plot(res_times, cstr_selectivity, 'g--', label='CSTR Selectivity')
        plt.title(f'Selectivity Comparison: {reaction_id} @ {temp}°C, C_A0={conc}')
        plt.xlabel('Residence Time (τ) (hours)'); plt.ylabel('Selectivity (B/I)')
        plt.legend(); plt.grid(True, linestyle='--'); plt.ylim(bottom=0)
        plt.savefig(os.path.join(output_folder, f'{reaction_id}_T{temp}_C{conc}_selectivity_comparison.png'))
        plt.close()

    pd.DataFrame(comparison_results).to_csv(FINAL_COMPARISON_SUMMARY, index=False, float_format='%.2f')
    print(f"\n\n{'='*50}\n✅ REACTOR COMPARISON COMPLETE\n{'='*50}")
    print(f"📊 A final summary with reactor recommendations has been saved to '{FINAL_COMPARISON_SUMMARY}'")

if __name__ == '__main__':
    compare_reactors()
