def l_Point(l,E0,p):
    cofactor=(p+1)/l
    while True:
        Q = E0.random_point()
        P = cofactor * Q                
        if P.order() == l:
            break
    return P
def walk(s_key,E,p):
    assert [e%(p+1)==0 for e in s_key]
    L=[l for l in prime_range(3,p) if (p+1)%l==0]
    for l,e in zip(L,s_key):
        for _ in range(e):
            P=l_Point(l,E,p)
            E=E.isogeny(P).codomain()
    #a2, a4 = E.a2(), E.a4()
    return E
s_key=[2,3,4];E=EllipticCurve(GF(419), [0,0,0,1,0]);E1=walk(s_key,E,419);M=E1.montgomery_model();print(M.a2()/M.a4())
