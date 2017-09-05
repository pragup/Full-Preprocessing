from help import Set_Covering_Algorithm as SC

preprocess_selection = 1
cleanup_selection = 1
preprocess = 1
SC("cap360.dat", preprocess, preprocess_selection, cleanup_selection, "greedy_score", khard=2)


