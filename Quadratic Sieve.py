from math import fabs,ceil, sqrt, exp, log
import random
from Factor import brent
from MillerRabin import is_probable_prime
from operator import add
from itertools import chain

def gcd(a,b): # Euclid's algorithm
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)


def prime_gen(n): # sieve of Eratosthenes, generates primes up to a bound n
    if n < 2:
        return []
    
    nums = []
    isPrime = []
    
    for i in range(0, n+1):#Creates list of numbers from 0 to n
        nums.append(i)
        isPrime.append(True)
        
    isPrime[0]=False
    isPrime[1]=False
    
    for j in range(2,int(n/2)):#tries all size gaps that make sense
        if isPrime[j] == True:
            for i in range(2*j,n+1,j):#starts from j+j, jumps by gap size j and crosses out that number
                isPrime[i] = False
                
    primes = []
    for i in range(0, n+1):#Adds leftovers
        if isPrime[i] == True:
            primes.append(nums[i])
            
    return primes


def pollard(N,factors):
    
    rem = N
    while True:
        if is_probable_prime(rem):
            factors.append(rem)
            break

        f = brent(rem)
        while f == rem:#ensures pollard rho returns a smaller factor
            f = brent(rem)
            
        if f and f < rem: #found a factor
            if is_probable_prime(f): #ensure f is prime
                #print("Pollard rho (Brent): Prime factor found: %s" % f)
                factors.append(f)
                rem = rem//f #other factor
            else: #factor is composite
                #print("Pollard rho (Brent): Non-prime factor found: %s" % f)
                rem_f = factor(f,factors) #recursive part
                rem = (rem//f) * rem_f #combines the two remainders
                factors.remove(rem_f)#removes tricky duplicate that got appended in 1st if stmt
        else: #no more factors found, rem is prime
            #print("No (more) small factors found.")
            break
                 
    return rem


def legendre(a, p): #legendre symbol of (a/p)
    return pow(a, (p - 1) // 2, p)
 
def tonelli(n, p): #tonelli-shanks to solve modular square root
    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        r = pow(n, (p + 1) // 4, p)
        return r,p-r
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return (r,p-r)


def smooth_bound(N): # finds optimal smoothness bound B
    B = exp(0.5*pow(log(N*log(log(N))),0.5))
    return int(B)

def find_base(N,P):
# generates a factor base with P primes
    global factor_base
    factor_base = []
    primes = prime_gen(P*10)# list to sieve for factors
    #print(primes)
    
    for p in primes:
    # find P primes such that N is a quadratic residue mod p
        if len(factor_base) == P:
            break
        if legendre(N,p) == 1:
            factor_base.append(p)
    return factor_base

def sieve_prep(N,interval=100):
# generates a sequence from Y(x) = x^2 - N, starting at x = ceil(sqrt(N)) 
    global sieve_seq
    sieve_seq = []
    for x in range(root,root+interval):
        sieve_seq.append(pow(x,2) - N)
    print(sieve_seq[18],sieve_seq[55],sieve_seq[83])
    return sieve_seq

def find_smooth(factor_base,N):
# tries to B-smooth numbers in sieve_seq, using sieving
    sieve_list = sieve_prep(N).copy() # the sieve list to be modified. We keep a copy of V for later
    
    print(len(sieve_list))
    print(str(sieve_list)+"\n")
    smooth_nums = []

    if 2 in factor_base:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i,len(sieve_list),2): # found the 1st even term, now every other term will also be even
            sieve_list[j] /= 2
        #print(sieve_list)
        
    for p in factor_base[1:]:
        residues = tonelli(N,p) #finds x such that x^2 = n (mod p). There are two start solutions except when p = 2
        print(residues)
        for r in residues: 
            for i in range((r-root)%p,len(sieve_list),p): # Now every pth term will also be divisible
                    sieve_list[i] = sieve_list[i]/p
                         
    print(sieve_list)
    print(sieve_list[18],sieve_list[55],sieve_list[83])
    xlist = [] #original x terms
    for i in range(len(sieve_list)):
        if sieve_list[i] == 1: # found B-smooth number
            smooth_nums.append(sieve_seq[i])
            xlist.append(i+root)
    return(smooth_nums,xlist)
                
    
def build_matrix(smooth_nums,factor_base):
# generates exponent vectors mod 2 from previously obtained smooth numbers, builds matrix
    matrix = []

    def factor(n,factor_base):#trial division from factor base
        factors = []
        for p in factor_base:
            if n % p == 0:
                factors.append(p)
        return factors

    for n in smooth_nums:
        exp_vector = [0]*len(factor_base)
        n_factors = factor(n,factor_base)
        print(n,n_factors)
        for i in range(len(factor_base)):
            if factor_base[i] in n_factors:
                exp_vector[i] = (exp_vector[i]+1)%2

        #print(n_factors, exp_vector)
        matrix.append(exp_vector)
    print("Matrix built:" + str(matrix))

    return(transpose(matrix))

    
def transpose(matrix):
#transpose matrix so columns become rows, makes list comp easier to work with
    new_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        new_matrix.append(new_row)
    return(new_matrix)

def gauss_elim(M):
#reduced form of gaussian elimination, finds rref and reads off the nullspace
#https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
    marks = [False]*len(M[0])
    
    for i in range(len(M)): #do for all rows
        row = M[i]
        print(row)
        for num in row: #search for pivot
            if num == 1:
                print("found pivot at column " + str(row.index(num)+1))
                j = row.index(num) # column index
                marks[j] = True 
                for k in chain(range(0,i),range(i+1,len(M))): #search for other 1s in the same column
                    if M[k][j] == 1:
                        for i in range(len(M[k])): M[k][i] = (M[k][i] + row[i])%2
                        print(M)
                break
    print(marks)
    
    #find dependent vectors to combine    
    M = transpose(M)
    print(M)
    free_rows = []
    for i in marks: #find free columns (which have now become rows)
        if i == False:
            free_rows.append(M[marks.index(i)])
    print(free_rows)

    solution_vec = []
    indices = []
    free_row = free_rows[0]
    for i in range(len(free_row)):
        if free_row[i] == 1: #rows with 1 in the same column will be dependent
            indices.append(i)
    for row in M:
        for i in indices:
            if row[i] == 1:
                solution_vec.append(M.index(row))
                break
            
    return(solution_vec)
    
def solve(solution_vec,smooth_nums,xlist,N):
    
    solution_nums = [smooth_nums[i] for i in solution_vec]
    print(solution_nums)

    Asquare = 1
    for n in solution_nums:
        Asquare *= n

    Bsquare = 1
    for n in xlist:
        Bsquare *= n**2

    print(str(Asquare)+" = "+str(Bsquare)+" mod "+str(N))
    a = sqrt(Asquare)
    b = sqrt(Bsquare)
    print(a,b)
    
    factor = gcd(fabs(b-a),N)
    
    print("Found factors! " + str(int(factor)) + ", " + str(int(N/factor)))

def QS(n,P):
#single polynomial version of quadratic sieve, up to the first P primes
    
    global N
    global root
    N = n
    root = ceil(sqrt(N))
    
    print(root)
    print("Attempting to factor {}...".format(N))
    
    print("Generating factor base...")
    factor_base = find_base(N,P) #generates a B-smooth factor base
    print(factor_base)

    B = factor_base[-1]
    print("Looking for {} {}-smooth relations...".format(len(factor_base)+1,B))
    smooth_nums,xlist = find_smooth(factor_base, N) #finds B-smooth relations, using sieving and Tonelli-Shanks
    print("Found {} smooth numbers.".format(len(smooth_nums)))
    print(smooth_nums)
    if len(smooth_nums) == 0:
        print("Not enough smooth numbers. Increase the sieve interval or size of the factor base.")
        quit()
    
    print("Building exponent matrix...")
    matrix = build_matrix(smooth_nums,factor_base) #builds exponent matrix mod 2 from relations

    print("Performing Gaussian Elimination...")
    solution_vec = gauss_elim(matrix) #solves the matrix for the null space, finds perfect squares
    print("Solution vector found: " + str(solution_vec))
    print("Solving congruence of squares...")
    solve(solution_vec,smooth_nums,xlist,N) #solves the congruence of squares to obtain factors

  
