# PerceivedControlIntegration


Figures 1A and 3A depict scripts titled "Rating_SERE.m" and "Rating_to_Play.m," respectively. These scripts are responsible for categorizing the confidence ratings provided by participants in the Rating task and the Play task, respectively.

Model:
"run_mle_uncertainty_condition.m": This script utilizes maximum likelihood estimation to fit four types of models, including the SE model, RE model, Integrated model, and Expected Value model. It incorporates "mle_modelFitting_simple.m" for this purpose.

Model Comparison:
"group_mle_condition_BIC.m": This script conducts model comparison among the four aforementioned models using the Bayesian Information Criterion (BIC).
"group_mle_condition.m": This script examines how participants weigh the SE and RE components under the SERE model.


Additionally, there are supplementary scripts:
"Rating_EvaGame.m": This script categorizes the likely ratings provided by participants in the Evaluation task.
