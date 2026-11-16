import json,os,urllib.request as u
E=os.environ.get('MAKEFUN_ENV_URL')
V=[0,100,320,330,500,900,0]
M=[0,0,8,6,2,3,-6]
J=[-21,-19,-12,-8,8,12,19,21]
I=[-11,-10,-9,-1,1,9,10,11]
G=[-11,-9,9,11]
O=[-10,-1,1,10]
SQ=[x for x in range(21,99)if x%10%9]
CT=[0]*120
for s in SQ:CT[s]=min(s%10-1,8-s%10)+min(s//10-2,9-s//10)
FLIP=bytes(i^8*(0<i<7 or 8<i<15)for i in range(256))
REV=[0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15]
CK={21:14,28:13,91:11,98:7}
TT={}
N=0;NB=2e5;BM=0;Z=119
class X(Exception):pass
def rot(B):return B[::-1].translate(FLIP)
def att(B,s):
 if B[s+9]==9 or B[s+11]==9:return 1
 for L,v in((J,10),(I,14)):
  for t in L:
   if B[s+t]==v:return 1
 for d in I:
  o=s+d
  while B[o]==0:o+=d
  if B[o]>6 and(B[o]&7)in(5,3+(d in O)):return 1
 return 0
def parse(f):
 a=f.split()
 W=bytearray([7]*120)
 for s in SQ:W[s]=0
 r=9;c=1;K=0;k=0
 for ch in a[0]:
  if ch=='/':r-=1;c=1
  elif ch<'A':c+=int(ch)
  else:
   q=r*10+c;v='PNBRQKpnbrqk'.index(ch)+1+2*(ch>'Z');W[q]=v
   if v==6:K=q
   elif v==14:k=q
   c+=1
 C=sum({'K':2,'Q':1,'k':8,'q':4}.get(ch,0)for ch in a[2])
 E=0 if a[3]=='-'else ord(a[3][0])+10*int(a[3][1])-86
 B=bytes(W)
 if a[1]=='b':
  B=rot(B);C=REV[C];E=Z-E;K,k=k and Z-k,K and Z-K
 return(B,C,E,K,k)
def mv(P,m):
 B,C,E,K,k=P
 f=m&127;t=(m>>7)&127;pr=(m>>14)&7
 p=B[f];q=B[t]
 W=bytearray(B)
 if p==1 and t==E and not q:W[t-10]=0
 W[t]=pr or p;W[f]=0
 if p==6:
  K=t
  if abs(t-f)==2:
   d=t>f;r=21+7*d;W[f+2*d-1],W[r]=W[r],0
 C&=CK.get(f,15)&CK.get(t,15)&(15-3*(p==6))
 E=(f+t)//2 if p==1 and abs(t-f)==20 else 0
 W=bytes(W)
 if att(W,K):return None
 return(rot(W),REV[C],Z-E,k and Z-k,K and Z-K)
def gen(P,cp,tt=0):
 B,C,E,K,k=P
 ms=[];A=ms.append
 for f in SQ:
  p=B[f]
  if not p or p>6:continue
  if p==1:
   o=f+10
   if B[o]==0:
    if f>80:
     for pr in(5,4,3,2):A((9e3,f|o<<7|pr<<14))
    elif not cp:
     A((0,f|o<<7))
     if f<39 and B[f+20]==0:A((0,f|(f+20)<<7))
   for o in(f+9,f+11):
    q=B[o]
    if q>7:
     if f>80:
      for pr in(5,4,3,2):A((9e3+V[q&7],f|o<<7|pr<<14))
     else:A((1e3+V[q&7],f|o<<7))
    elif o==E:A((1e3,f|o<<7))
  else:
   for d in(J,G,O,I,I)[p-2]:
    o=f+d
    while B[o]==0:
     if not cp:A((0,f|o<<7))
     if p%4==2:break
     o+=d
    if B[o]>7:A((1e3+V[B[o]&7],f|o<<7))
 if not cp and not att(B,K):
  k=K
  for s,R,b,x in((1,28,2,27),(-1,21,1,22)):
   if C&b and B[k+s]==B[k+2*s]==B[x]==0 and B[R]==4 and att(B,k+s)==att(B,k+2*s)==0:A((0,k|(k+2*s)<<7))
 ms.sort(key=lambda x:x[0]+(x[1]==tt)*1e5,reverse=1)
 return[m for _,m in ms]
def ev(B):
 s=0
 for q in SQ:
  p=B[q]
  if p:
   t=p&7;v=V[t]+M[t]*CT[q]
   s+=v if p<8 else-v
 return s
def qs(P,a,b,pl):
 global N
 N+=1
 if N>NB:raise X
 B,C,E,K,k=P
 inchk=att(B,K)
 st=ev(B)
 if not inchk:
  if st>=b:return b
  a=max(a,st)
 n=0
 for m in gen(P,0 if inchk else 1):
  Q=mv(P,m)
  if not Q:continue
  n+=1
  s=-qs(Q,-b,-a,pl+1)
  if s>=b:return b
  a=max(a,s)
 if inchk and not n:return-1e5+pl
 return a
def sr(P,d,a,b,pl):
 global N,BM
 N+=1
 if N>NB:raise X
 B,C,E,K,k=P
 if d<=0:return qs(P,a,b,pl)
 key=(B,C,E)
 if pl and key in H:return 0
 e=TT.get(key);tt=0
 if e and e[0]>=d and e[2]==0:return e[1]
 if e:tt=e[3]
 bm=0;n=0
 for m in gen(P,0,tt):
  Q=mv(P,m)
  if not Q:continue
  n+=1
  s=-sr(Q,d-1,-b,-a,pl+1)
  if s>a:a=s;bm=m;BM=m*(not pl)or BM
  if s>=b:
   TT[key]=(d,1,s,m)
   return b
 if not n:return-1e5+pl if att(B,K)else 0
 TT[key]=(d,1-(bm!=0),a,bm)
 return a
def post(p,b):
 return json.load(u.urlopen(E+p,json.dumps(b).encode()))
def choose(g):
 P=parse(g['fen'])
 L=g['legal_moves']
 H.add(P[:3])
 global N,BM
 N=0;BM=0
 for d in range(1,99):
  try:sr(P,d,-1e6,1e6,0)
  except X:break
 if not BM:return L[0]
 f=BM&127;t=(BM>>7)&127;pr=BM>>14
 if ' b ' in g['fen']:f,t=Z-f,Z-t
 s=''.join(chr(96+v%10)+str(v//10-1)for v in(f,t))
 if pr:s+='nbrq'[pr-2]
 if s not in L:return L[0]
 return s
while True:
 g=post('/game',{});H=set()
 if g['status']=='done':break
 while g['status']=='your_move':
  F=g['legal_moves'][0]
  try:m=choose(g)
  except:m=F
  for v in(m,F):
   try:g=post('/move',{'game':g['game'],'move':v});break
   except:pass
  else:exit(1)
