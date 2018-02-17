import numpy as np
import ast

def cliquefinder(i,B,v):
    x = B[i]==1
    M = B[x,:]
    M = M[:,x]
    for j in range(v):
        x = M[j]==1
        M = M[x,:]
        M = M[:,x]
        if M.shape[0]==j+1:
            return M

def buildadj(file,l):
    f = open(file,'r')
    S = [set(ast.literal_eval(line)) for line in f]
    f.close()
    A = np.zeros((len(S),len(S)))
    for i,entry1 in enumerate(S):
        for j,entry2 in enumerate(S):
            A[i,j] = len(entry1.intersection(entry2))
    A[A!=l] = 0
    A[A==l] = 1
    np.save(file[:-4]+'graph.npy',A)
    return A
        
def forcliquer(A,outfile):
    f = open("../cliquer-1.21/{}".format(outfile),'w')
    verts = A.shape[0]
    edges = int(sum(sum(A)))
    f.write("p edge {} {}\n".format(verts,edges))
    for i in range(verts):
        for j in range(verts):
            if i<j:
                if A[i,j] == 1:
                    f.write("e {} {}\n".format(i+1,j+1))

def shift(a,b,ring):
    c = []
    for i,char in enumerate(ring):
        c.append((a[i]+b[i])%int(char))
    return tuple(c)

def shiftcheck(ring):
    if ring == '66':
        Pts = {(i,j) for i in range(6) for j in range(6)};
        f = open("z6z6.txt",'r')
    elif ring == '433':
        Pts = {(i,j,k) for i in range(4) for j in range(3) for k in range(3)};
        f = open("z4z3z3.txt",'r')
    else:
        print("The specified ring is not yet implemented.")
    blocks = [set(ast.literal_eval(line)) for line in f]
    wins = [1 for i in range(len(blocks))]
    for j,block in enumerate(blocks):
        listks = ['' for i in range(36)]
        for k,point in enumerate(Pts):
            newblk = {shift(blkpt,point,ring) for blkpt in block}
            listks[k] = len(newblk.intersection(block))
            if listks[k]==6:
                wins[j]+=1
    return wins

def spence(P,pshift,K,ktype):
    H1 = {(0,0),(1,0),(2,0)}
    H2 = {(0,0),(0,1),(0,2)}
    H3 = {(0,0),(1,1),(2,2)}
    H4 = {(0,0),(1,2),(2,1)}
    Pts = H1.union(H2).union(H3).union(H4)
    planes = [H1,H2,H3,H4]
    k = [0,1,2,3]
    totalPts = sorted([(h[0],h[1],k[i]) for h in Pts for i in range(4)])
    D = {(shift(h,pshift[i],'33')[0],shift(h,pshift[i],'33')[1],k[K[i]]) for i in range(1,4) for h in planes[P[i]]}.union({(h[0],h[1],k[K[0]]) for h in Pts.difference(planes[P[0]])})
    count = [0 for i in range(36)]
    for d1 in D:
        for d2 in D:
            count[totalPts.index(spencediff(d1,d2,ktype))]+=1
    return count,sorted(list(D))
            

def spencediff(a,b,ktype):
    if ktype=='z4':
        d = ((a[0]-b[0])%3,(a[1]-b[1])%3,(a[2]-b[2])%4)
    elif ktype == 'klein':
        klein = [(0,0),(0,1),(1,0),(1,1)]
        x = klein[a[2]]
        y = klein[b[2]]
        c = ((x[0]+y[0])%2,(x[1]+y[1])%2)
        d = ((a[0]-b[0])%3,(a[1]-b[1])%3,klein.index(c))
    return d
