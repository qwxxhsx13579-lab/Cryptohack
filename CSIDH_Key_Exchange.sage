ells = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 
        71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 
        149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 
        227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
        307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587]
p = 4 * prod(ells) - 1
F = GF(p)
E0 = EllipticCurve(F, [1, 0])
a_priv = [-1, -2, -3, -3, -2, -3, -3, 0, 2, -1, 2, -1, -2, -3, 1, 2, 1, 2, 0, 0, 1, -1, 0, 2, -1, 0, 0, 0, 1, -1, -3, 1, -1, -3, -3, 2, 2, 1, -1, -1, 1, 0, 1, 1, 1, -2, 2, 2, -2, -2, 0, 0, 2, 0, -1, -3, -2, -2, 0, -1, -3, -1, -2, -3, -2, 2, 1, 1, -2, 0, 1, -1, -3, 2]
b_priv = [-1, -1, 0, 1, 2, 0, 2, -1, -3, 1, 0, -2, -2, 2, -1, -2, -3, -3, -3, 2, 2, 2, -2, -1, 1, -2, 0, -3, -1, 1, -1, -1, -3, -1, -2, 1, -1, -2, -3, 1, 0, -1, 1, 2, 2, 0, 0, -1, -2, -2, 1, -1, 1, 1, 1, 1, 0, 0, 0, -3, -2, -1, 2, 0, -3, -2, 1, 1, -2, -1, -1, 2, 0, 1]
def l_Point(l,E,p):
    cofactor=(p+1)/l
    while True:
        Q = E.random_point()
        P = cofactor * Q                
        if P != E(0):
            break
    return P
def walk(s_key,E,p):
    #assert [e%(p+1)==0 for e in s_key]
    L=ells;assert len(s_key) == len(ells)
    for l,e in zip(L,s_key):
        if e>0:
            for _ in range(abs(e)):
                P=l_Point(l,E,p)
                E=E.isogeny(P).codomain()
        elif e<0:
            Et = E.quadratic_twist(2)
            for _ in range(-e):
                P=l_Point(l,Et,p)
                Et=Et.isogeny(P).codomain()
        E = Et.quadratic_twist(2)
     
    return E             
print("hhh")
Ea=walk(a_priv,E0,p)#print(Ea)
Eb=walk(b_priv,E0,p)
print(walk(b_priv, Ea, p).is_isomorphic(walk(a_priv, Eb, p)))
Eab=walk(b_priv,Ea,p);M=Eab.montgomery_model();print(M.a2()/M.a4())
