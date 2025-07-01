
### This is a bare bones implementation of CCA + permutation testing. Please be aware that your X and Y data variables may not start at 8.
## You may also want a different number of permutations and components to evaluate 


# Canonical Correlation Analysis with Permutation Testing
n_permutations = 1000
n_components = 3

# Select variables for X and Y datasets
X = dfx.iloc[:, 8:]
Y = dfy.iloc[:, 8:]

# Create array for permutation values
permutation_scores = np.zeros((n_components, n_permutations))

# Begin CCA and Permutation Testing
for i in range(n_components):
    est.fit(X, Y)
    X_orig, Y_orig = est.transform(X, Y)
    correlation = np.corrcoef(X_orig[:, i], Y_orig[:, i])[0, 1]

    for j in range(n_permutations):
        #Only Y dataset is shuffled
        Y_shuffle = shuffle(Y, random_state = j)
        X_perm, Y_perm = est.transform(X, Y_shuffle)
        permutation_scores[i, j] = np.corrcoef(X_perm[:, i], Y_perm[:, i])[0, 1]
    # Compare true correlation to correlation derived from shuffling data
    p_value = (np.sum(permutation_scores[i] >= correlation) + 1)/(n_permutations +1)

    print(f"Component {i + 1}:")
    print(f"Observed Correlation: {correlation:.3f}")
    print(f"Permuted Correlation Range: {permutation_scores.min():.3f} to {permutation_scores.max():.3f}")
    print(f"P-value: {p_value:.3f}")
