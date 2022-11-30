import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest
import argparse
from statsmodels.stats.multitest import multipletests


def check_intervals_intersect(first_ci, second_ci):   
    return second_ci[0] < first_ci[0] < second_ci[1] or second_ci[0] < first_ci[1] < second_ci[1] or first_ci[0] < second_ci[0] < first_ci[1] or first_ci[0] < second_ci[1] < first_ci[1]


def check_dge_with_ci(first_table, second_table):
    ci_test_results = []
    for name in first_table.columns:
      first_CI = st.t.interval(alpha=0.95, # 95% доверительный интервал
              df=len(first_table[name]) - 1, # число степеней свободы - 1
              loc=np.mean(first_table[name]), # Среднее
              scale=st.sem(first_table[name]))
      second_CI = st.t.interval(alpha=0.95, # 95% доверительный интервал
              df=len(second_table[name]) - 1, # число степеней свободы - 1
              loc=np.mean(second_table[name]), # Среднее
              scale=st.sem(second_table[name]))
      ci_test_results.append(not(check_intervals_intersect(first_CI, second_CI)))
    return ci_test_results


def check_dge_with_ztest(first_table, second_table):
    z_test_results = []
    for name in first_table.columns:
      z_test = ztest(first_table[name],second_table[name])
      z_test_list = [z_test[1] < 0.05, z_test[1]]
      z_test_results.append(z_test_list)
    return z_test_results


def check_dge_with_ztest(first_table, second_table):
    z_test = ztest(first_table, second_table)
    return z_test

    
def check_dge_with_ttest(first_table, second_table):
    t_test = st.ttest_ind(first_table, second_table)
    return t_test


def dge(first_cell_type_expressions_path, second_cell_type_expressions_path, save_results_table, test_method='t', p_cor_method='hs', alpha=0.05):
    expression_data_1 = pd.read_csv(first_cell_type_expressions_path, index_col=0)
    expression_data_2 = pd.read_csv(second_cell_type_expressions_path, index_col=0)
    
    mean_diff = np.mean(expression_data_1) - np. mean(expression_data_2)
    
    ci_test_results = check_dge_with_ci(expression_data_1, expression_data_2)

    if test_method == 'z':
        test_results = check_dge_with_ztest(expression_data_1, expression_data_2)
    else:
        if test_method != 't':
            print("I don't know this test, so I use t-test")
            test_method = 't'
        test_results = check_dge_with_ttest(expression_data_1, expression_data_2)
        
    p_adjusted = multipletests(test_results[1], method = p_cor_method)
    test_output = p_adjusted[1] < alpha
    results = {
        'ci_test_result' : ci_test_results,
        "mean_diff": mean_diff,
        f"{test_method}_test_p_values": test_results[1],
        f"{test_method}_test_p_adjusted": p_adjusted[1],
        f"{test_method}_test_results": test_output
        }
    results = pd.DataFrame(results)
    results.to_csv(f'{save_results_table}')

parser = argparse.ArgumentParser(description='Tool for differential gene expression')
parser.add_argument("--first_cell_type_expressions_path", type=str, help='Path to table with expressions of first cell type', required=True)
parser.add_argument('--second_cell_type_expressions_path', type=str, help='Path to table with expressions of second cell type', required=True)
parser.add_argument('--save_results_table', type=str, help='Output table name', required=True)
parser.add_argument('--test_method', type=str, help='DGE method testing. Choose "z" for z-test or "t" for t-test (default)', default='t')
parser.add_argument('--p_adjusted_method', type=str, help='Choose the method for muptiple hypothesis testing. "holm-sidak" by default. See more https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html', default='hs')
parser.add_argument('--alpha', type=float, help='Choose the level of DGE testing. alpha = 0.05 by default', default=0.05)

args = parser.parse_args()

dge(args.first_cell_type_expressions_path, args.second_cell_type_expressions_path, args.save_results_table, args.test_method, args.p_adjusted_method, args.alpha)

