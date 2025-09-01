# Reaction Kinetics Analysis  

This repository is dedicated to the **analysis of chemical reaction kinetics**.  
It provides a structured environment for:  

- Storing experimental data  
- Running Python scripts for data processing and calculations  
- Collecting and visualizing generated plots  

The project is organized by problem, allowing for independent analysis of different datasets.  

---

## 📂 Repository Structure  

The repository is organized into distinct directories for each analysis problem, along with a central location for key plots and scripts.  

### Main Directory  

- **Problem 1/** → Contains all the data, scripts, and plots related to the first analysis problem.  

---

## 📁 Contents of Problem 1  

The **Problem 1** directory is organized into subfolders and files corresponding to the kinetic analysis workflow.  

### 🔹 Data and Plots  

- **RX1/**, **RX2/**, **RX3/** → Subfolders containing plots and analysis outputs for each reaction.  
- **RX1.csv**, **RX2.csv**, **RX3.csv** → Raw experimental data files for each reaction.  
- **RX1_reaction_summary.csv**, **RX1_series_reaction_summary.csv** → Summary data files generated from the analysis of Reaction 1.  
- **all_reactions_summary.csv** → A consolidated summary of key findings from all reactions.  

### 🔹 Python Scripts  

- **arrhenius.py** → Calculations related to the Arrhenius equation.  
- **diff_rateanalyzer.py** → Differential rate analysis.  
- **integral_rateanalyzer.py** → Integral rate analysis.  
- **Rate_constant.py** → Rate constant calculations.  
- **RX_series_plot.py** → Plot generation for reaction series.  
- **reaction_order_summary.csv** → Summary of reaction order determination.  

---

## ⚙️ Usage  

To use the analysis scripts and regenerate plots:  

1. **Clone the repository**  
   ```bash
   git clone https://github.com/your-username/your-repository-name.git
   cd your-repository-name
   ```

2. **Install dependencies**  
   The scripts require Python and standard scientific libraries.  
   ```bash
   pip install pandas numpy matplotlib
   ```

3. **Run the analysis**  
   Example:  
   ```bash
   python arrhenius.py
   python RX_series_plot.py
   ```

---

