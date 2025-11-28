# Energy System Resilience Optimization Task

## 1. Introduction and Objective
Assess multi-bus network resilience under Dunkelflaute and thermal outage using LP.

## 2. Model Structure and Assumptions 💡
### Modeling Simplifications
* **Transport Model Used:** Simple nodal power balance is used, ignoring AC/DC power flow, losses, and voltage limits.
* **Linearity:** All costs and dynamics (SOC, efficiency) are assumed perfectly linear.
* **No Start-up Costs:** Thermal generation has a fixed marginal cost (£50/MWh), assuming zero start-up/shutdown costs.

### Critical Assumption: Perfect Foresight
* **The model is deterministic:** It solves the entire 672-hour period simultaneously and assumes **perfect foresight** of all load, wind, and outage events from $t=0$.

## 3. Scenario Definitions and Results
### Scenario Definitions

| Scenario | Name | Time Horizon | Load & Wind Profile | Thermal Capacity Constraint |
| :---: | :--- | :--- | :--- | :--- |
| **S1** | Normal Conditions | 672 Hours | Baseline ($100\%$ Load, $100\%$ Wind) | $\text{P}_{\text{Th}} \le 250\text{ MW}$ (Always Available) |
| **S2** | Dunkelflaute Baseline | 672 Hours | **Week 2:** Load $\times 1.40$, Wind $\times 0.40$ | $\text{P}_{\text{Th}} \le 250\text{ MW}$ (Always Available) |
| **S3** | Planned 24h Outage | 672 Hours | Same as S2 | $\text{P}_{\text{Th}} \le 0\text{ MW}$ for Hours 216-240 (Known) |
| **S4** | Unplanned 24h Outage | 672 Hours | Same as S2 | $\text{P}_{\text{Th}} \le 0\text{ MW}$ for Hours 216-240 (Treated identically to S3) |
Note on S3 and S4: In the deterministic LP model, S3 and S4 are mathematically equivalent because the optimizer has perfect foresight of the $\text{P}_{\text{Th,Max}}=0$ constraint from the start of the simulation.

## 4. Discussion: Why S2, S3, and S4 Results are Identical 📢
* **Methodological Reason (S3 vs. S4):** The model's deterministic nature means the **Planned Outage (S3)** and the **Unplanned Outage (S4)** are mathematically identical, as the optimizer has foreknowledge of the constraint in both cases.
* **Technical Reason (S2 vs. S3/S4):** The identical results suggest the $24\text{-hour}$ thermal outage constraint was **non-binding (redundant)**. The system was already so severely stressed by the Dunkelflaute in S2 that the maximum unserved load was limited by the **overall energy deficit**, not the short-term loss of the thermal plant.

## 5. Next Steps and Future Work
* **To capture true resilience:** The model should be re-run as a **Sequential Dispatch Model** or a **Multi-Stage Stochastic Program** to introduce uncertainty and prevent perfect foresight.
