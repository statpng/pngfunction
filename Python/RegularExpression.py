import re

dataset2 = pd.DataFrame( np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
              columns=['A', "T_1", "T_2"] )

# [ i if re.search("T_", i) else "X" for i in dataset.columns ]
[ i for i in dataset.columns if re.search("T_", i) ]
[ i for i in range(dataset2.shape[1]) if re.search("T_", dataset2.columns[i]) ]

