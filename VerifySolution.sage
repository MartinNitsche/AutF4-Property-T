# This SAGE-script verifies with interval arithmetic that up to a sufficiently small error
# 1/48*Laplacian^2 - eps*Laplacian = W := sum w_i^*w_i in R[SAut(F_4)],
# where eps > 0 and the w_i are read from the files "support.dat" and "coefficients.dat"
# that can be downloaded (in compressed form) from https://zenodo.org/record/7065231
# This proves that Aut(F_4) has property (T).
#
# The generators of F_4 are represented by the numbers 1-4, their inverses by the negatives.
# The generators of SAut(F_4) are represented by tuples (i,j) with 1 <= |i|, |j| <= 4.
# They act on F_4 by right-multiplying the i-th generator with the j-th generator of F_4.
# Elements of SAut(F_4) are represented either as words over its generators or as a string
# of a 4-tuple containing the value of the automorphism on each generator of F_4.
# We think of SAut(F_4) as acting on F_4 from the right.
#
# This script is not optimized for efficiency and may take many hours to run.

S = {} # generators of SAut(F_4), indexed by a string representation
for i in [-4,-3,-2,-1,1,2,3,4]:
  for j in [-4,-3,-2,-1,1,2,3,4]:
    if abs(i) != abs(j):
      S[str(i)+"/"+str(j)] = (i,j)

# takes two words a,b over S and returns a^{-1}b as a word over S
def differenceWord(a, b):
  return [(l[0],-l[1]) for l in reversed(a)] + b

# takes a,b in F_4 as reduced words and returns a*b as a reduced word
def multiplyInF4(a, b):
  result = list(a)
  for letter in b:
    if result and letter == -result[-1]: # reduce y*x*x^{-1} to y
      result.pop()
    else:
      result.append(letter)
  return tuple(result)

# takes a word in S and returns a string representation of the automorphism in 4-tuple form
def automorphism(word):
  tmp = [(1,), (2,), (3,), (4,)] # start with the identity automorphism
  for g in reversed(word): # inductively write result = shorterWord * tmp
    if g[0]*g[1] > 0:
      multiplyWith = tmp[abs(g[1])-1]
    else:
      multiplyWith = [-k for k in reversed(tmp[abs(g[1])-1])]
    if g[0] > 0:
      tmp[abs(g[0])-1] = multiplyInF4(tmp[abs(g[0])-1], multiplyWith)
    else:
      tmp[abs(g[0])-1] = multiplyInF4(multiplyWith, tmp[abs(g[0])-1])
  return str(tuple(tmp)) # finally, result = emptyWord * tmp = tmp

print("please run with >= 6 GB RAM. reading data...")
support = []  # support of the w_i, as words over S
mat_wi = None # matrix with rows corresponding to the elements w_i

with open("support.dat", "r") as file:
  for line in file:
    support.append([S[string] for string in line.rstrip(", \n").split(", ")])
support.append([]) # add neutral element at the very end
assert(max(len(word) for word in support) <= 3) # used in the error estimate at the end

with open("coefficients.dat", "r") as file:
  mat_wi = matrix(RIF, sum(1 for line in file), len(support)) # initialize to 0-matrix
with open("coefficients.dat", "r") as file:
  for linenumber,line in enumerate(file): # every line contains the coefficients of one w_i
    coefficients = [RealInterval(float(s)) for s in line.rstrip(", \n").split(",")]
    # last coefficient is for the neutral element, ensuring w_i is in augmentation ideal
    coefficients.append(-sum(coefficients))
    assert len(coefficients) == len(support)
    for i in range(len(support)):
      mat_wi[linenumber,i] = coefficients[i]

# estimate running time
import time
time0 = time.process_time()
_ = mat_wi[:400,:400].T*mat_wi[:400,:400] # only needed for run time estimate
time_est = (time.process_time() - time0) / 400^3 * mat_wi.ncols()^2 * mat_wi.nrows() * 1.1
print("main computation... (estimated running time: %d cpu minutes)" % round(time_est/60))

mat_W = mat_wi.T*mat_wi # W will be the group ring element associated to this positive matrix
print("finishing up...")

W = {} # coefficients of W as intervals, indexed by group elements in 4-tuple form
for i in range(len(support)):
  for j in range(len(support)):
    difference = automorphism(differenceWord(support[i], support[j]))
    coefficient = W.get(difference, RealInterval(0)) # if no entry yet, start with 0
    W[difference] = coefficient + mat_W[i,j]

# from W now subtract 1/48*Laplacian^2 - eps*Laplacian, then compare to 0
laplacian = {(s,):RealInterval(-1) for s in S.values()}
laplacian[()] = -sum(laplacian.values()) # neutral element coefficient
for leftWord in laplacian:
  for rightWord in laplacian:
    difference = automorphism(differenceWord(list(leftWord), list(rightWord)))
    assert difference in W
    W[difference] -= laplacian[leftWord] * laplacian[rightWord] / RealInterval(48)
eps = W[automorphism([S["1/2"]])]
for word in laplacian:
  W[automorphism(word)] += laplacian[word] * eps
xNorm = RealInterval(0) # 1-norm of the error
for coefficient in W.values():
  xNorm += abs(coefficient)

# check that the error satisfies the estimate of Netzer-Thom, Lemma 2.1 with d=3
print("1/48*Laplacian^2 - eps*Laplacian >= 0 for eps >= %s" % (eps-16*xNorm))
print("success!" if eps > 16*xNorm else "verification failed")
