import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- Configuration ---
REACTION_FILES = ['RX1.csv', 'RX2.csv', 'RX3.csv']
MASTER_SUMMARY_FILE = 'reaction_type_summary.csv'
# Threshold: If the final concentration of B is less than 95% of its peak
# concentration, we classify it as a series reaction.
SERIES_REACTION_THRESHOLD = 0.95

def auto_classify_reaction_type():
    """
    Analyzes a list of reaction data files to automatically determine if each
    reaction is series or parallel based on the concentration profile of product B.
    Generates plots and a summary file with the classification.
    """
    print("🚀 Starting Automatic Reaction Type Classifier...")
    
    master_summary_list = []

    # --- Main loop to process each CSV file ---
    for csv_file in REACTION_FILES:
        reaction_id = os.path.splitext(csv_file)[0]
        output_folder = f'{reaction_id}_type_classification_plots'
        os.makedirs(output_folder, exist_ok=True)
        
        print(f"\n{'='*50}\nProcessing Reaction File: {csv_file}\n{'='*50}")

        try:
            df = pd.read_csv(csv_file)
        except FileNotFoundError:
            print(f"❌ WARNING: The file '{csv_file}' was not found. Skipping.")
            continue
        
        # Loop through each individual experiment in the file
        for temp in df['Temperature'].unique():
            print(f"  - Analyzing Temperature: {temp}°C")
            temp_df = df[df['Temperature'] == temp]
            
            for conc in temp_df['Initial Concentration'].unique():
                subset_df = temp_df[temp_df['Initial Concentration'] == conc].copy()
                
                # --- The Core Classification Logic ---
                conc_b_max = subset_df['B'].max()
                conc_b_final = subset_df['B'].iloc[-1]
                
                # Avoid division by zero if no product is formed
                if conc_b_max > 0 and (conc_b_final / conc_b_max) < SERIES_REACTION_THRESHOLD:
                    reaction_type = "Series"
                    explanation = "([B] decreases after peak)"
                else:
                    reaction_type = "Parallel"
                    explanation = "([B] rises and plateaus)"
                
                print(f"    - C_A0={conc}: Determined Type = {reaction_type} {explanation}")

                # Append the result to our master list
                master_summary_list.append({
                    'Reaction_ID': reaction_id,
                    'Temperature_C': temp,
                    'Initial_Conc_A_mg/ml': conc,
                    'Determined_Reaction_Type': reaction_type,
                    'Peak_Conc_B': conc_b_max,
                    'Final_Conc_B': conc_b_final
                })

                # --- Plotting for Visual Confirmation ---
                plt.figure(figsize=(12, 8))
                plt.plot(subset_df['Time'], subset_df['A'], 'o-', color='blue', label='Reactant [A]')
                plt.plot(subset_df['Time'], subset_df['B'], 's-', color='green', label='Product [B]')
                plt.plot(subset_df['Time'], subset_df['I'], '^-', color='red', label='Impurity [I]')
                
                # Highlight the max and final points of B
                plt.plot(subset_df['Time'].iloc[-1], conc_b_final, 'x', color='darkgreen', markersize=12, label=f'Final [B] = {conc_b_final:.2f}')
                plt.plot(subset_df.loc[subset_df['B'].idxmax(), 'Time'], conc_b_max, '*', markersize=18, color='gold', 
                         markeredgecolor='black', label=f'Max [B] = {conc_b_max:.2f}')
                
                plt.title(f'Determined Reaction Type: {reaction_type}\n{reaction_id} @ {temp}°C (Initial [A] = {conc})', fontsize=16)
                plt.xlabel('Time (h)')
                plt.ylabel('Concentration (mg/ml)')
                plt.legend()
                plt.grid(True, linestyle='--')

                filename = f'{reaction_id}_T{temp}_C{conc}_classification.png'
                filepath = os.path.join(output_folder, filename)
                plt.savefig(filepath)
                plt.close()

    # --- Save the final consolidated summary ---
    if master_summary_list:
        results_df = pd.DataFrame(master_summary_list)
        results_df.to_csv(MASTER_SUMMARY_FILE, index=False, float_format='%.4f')
        print(f"\n\n{'='*50}\n✅ BATCH CLASSIFICATION COMPLETE\n{'='*50}")
        print(f"📊 A consolidated summary of reaction types has been saved to '{MASTER_SUMMARY_FILE}'")
        print(f"📈 All classification plots have been saved in their respective folders.")
    else:
        print("\n\nNo data was processed. Please check if your CSV files exist and are named correctly.")

if __name__ == '__main__':
    auto_classify_reaction_type()
