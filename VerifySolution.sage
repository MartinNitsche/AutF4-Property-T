# This SAGE-script verifies with interval arithmetic that up to a sufficiently small error
# 1/48*Laplacian^2 - eps*Laplacian = W := sum w_i^*w_i in R[SAut(F_4)],
# where eps > 0 and the w_i are read from the files "support.dat" and "coefficients.dat"
# that can be downloaded (in compressed form) from https://zenodo.org/record/7065231
# This proves that Aut(F_4) has property (T).
#
# The generators of F_4 are represented by the numbers 1-4, their inverses by the negatives.
# The generators of SAut(F_4) are represented by tuples (i,j) with 1 <= |i|, |j| <= 4.
# They act on F_4 by right-multiplying the i-th generator with the j-th generator of F_4.
# Elements of SAut(F_4) are represented either as words over its generators
# or as a 4-tuple containing the value of the automorphism on each generator of F_4.
# We think of SAut(F_4) as acting on F_4 from the _right_.
# The computation may take many hours and may require up to 8 GB of RAM plus swap space.

S = {} # generators of SAut(F_4), indexed by a string representation
for i in [-4,-3,-2,-1,1,2,3,4]:
  for j in [-4,-3,-2,-1,1,2,3,4]:
    if abs(i) != abs(j):
      S[str(i)+"/"+str(j)] = (i,j)

support = [] # support of the w_i, as words over S
with open("support.dat", "r") as file:
  for line in file.readlines():
    support.append([S[string] for string in line.rstrip(", \n").split(", ")])
support.append([]) # add neutral element
assert(max(len(word) for word in support) <= 3) # used in the error estimate at the end

summands = [] # coefficients of the w_i, approximated by real intervals
with open("coefficients.dat", "r") as file:
  for line in file.readlines():
    coefficients = [RealInterval(float(s)) for s in line.rstrip(", \n").split(",")]
    coefficients.append(-sum(coefficients)) # for neutral element
    assert len(coefficients) == len(support)
    summands.append(coefficients)

def differenceWord(a, b): # returns a^{-1}b as a word over S
  return [(l[0],-l[1]) for l in reversed(a)] + b

def multiplyInF4(a, b): # returns a*b as a reduced word over the generators of F_4
  result = list(a)
  for letter in b:
    if result and letter == -result[-1]: # reduce y*x*x^{-1} to y
      result.pop()
    else:
      result.append(letter)
  return tuple(result)

def automorphism(word): # returns the automorphism in 4-tuple form
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
  return tuple(tmp) # finally result = emptyWord * tmp = tmp

W = {} # coefficients of W as intervals, indexed by group elements in 4-tuple form
for i in range(len(support)): # the order of the 3 for-loops is switched for efficiency
  print("progress: %d / %d" % (i+1, len(support)))
  for j in range(len(support)):
    difference = automorphism(differenceWord(support[i], support[j]))
    coefficient = W.get(difference, RealInterval(0)) # if no entry yet, start with 0
    for summand in summands:
      coefficient += summand[i] * summand[j]
    W[difference] = coefficient

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
xNorm = 0 # 1-norm of the error
for coefficient in W.values():
  xNorm += abs(coefficient)

# check that the error satisfies the estimate of Netzer-Thom, Lemma 2.1 with d=3
print("1/48*Laplacian^2 - eps*Laplacian >= 0 for eps >= %s" % (eps-16*xNorm))
print("success!" if eps > 16*xNorm else "verification failed")
