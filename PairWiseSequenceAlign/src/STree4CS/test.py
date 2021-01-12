from STree4CS import STree4CS as STree

# Suffix-Tree example.
st = STree.STree4CS([[1, 10, 5, 3, 200, 8, 10, 2]])
print(st.find([1, 10])) # 0
print(st.find_all([10])) # [1, 6]

# Building a generalized Suffix-Tree.
A = [[1, 2, 3], [4, 5, 6, 2, 3, 7], [1, 2, 3, 4]]
st = STree.STree4CS(A)

# print the longest common subsequence for the set A
print(st.lcs()) # [2, 3]

# Sequence Covering similarity example
S=[[1,1,2,2,3,4,1,1,5,6], [1,2,4,3,4,5,7,5,1], [6,5,1,7,4,5,6]]
s=[1,1,5,7,5,1,7,4]

# Build the generalized Suffix-Tree for S
st = STree.STree4CS(S)

# get the S-optimal covering similarity for s
score, lbreak, lss = st.evaluateDichotomic(s)
print(score) # 0.75
print(lss) # the optimal covering [[1, 1, 5], [7, 5, 1], [7, 4]], 