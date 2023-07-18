- ```python
  #Solve the equilibrium problem for a heterodimer, two homodimers, and two unimolcular structures
  #Solve the equilibrium iteratively.
  #This is, so far, bare bones.  The user needs to enter ct (total strand concentrations) and
  #the 5 folding free energy changes
  #Different starting points are tried, but the user must look at the results to see what solutions
  #, if any, are physical.
  
  from scipy.optimize import fsolve
  import math
  
  # assume experiment at 37 degrees=310.15K
  T = 310.15
  # gas constant in kcal
  R = 1.9872036e-3
  
  # constants to be set, enter -DG in kcal/mol
  # ct is the total concentrations of strands in M.  Right now A and B are assumed to be mixed 1:1
  ct = 1e-6
  # Here are the equilibrium constants.  For now, I have entered free energy changes in kcal/mol.
  # These expressions convert the free energy change to equilibrium constant.
  KAB = math.exp(21 / (R * T))
  KA2 = math.exp(14.4 / (R * T))
  KB2 = math.exp(17.1 / (R * T))
  KA = math.exp(8 / (R * T))
  KB = math.exp(19.2 / (R * T))
  
  print("KAB = "+str(KAB))
  print("KA2 = "+str(KA2))
  print("KB2 = "+str(KB2))
  print("KA = "+str(KA))
  print("KB = "+str(KB))
  
  def func(x):
  
  
      #species
      #x[0] = [AB]
      #x[1] = [A]
      #x[2] = [B]
      #x[3] = [A2]
      #x[4] = [B2]
      #x[5] = [A']
      #x[6] = [B']
  
  </mark>
      return[x[0]/(x[1]*x[2])-KAB,x[3]/(x[1]*x[1])-KA2,x[4]/(x[2]*x[2])-KB2,x[5]/x[1]-KA,x[6]/x[2]-KB,\
             x[1]+x[5]+2*x[3]+x[0]-ct/2,x[2]+x[6]+2*x[4]+x[0]-ct/2]
  </mark>  
    
    
  #good start for good config [1e-5,1e-10,1e-10,1e-8,1e-8,1e-8,1e-8]
  
  #Solve using Newton's method.  This tends to get stuck in locally optimal solutions.
  #In particular, a solution with a negative root for any of the concentrations is not physical.
  #Or a solution with the sum of concentrations>Ct is also not physical.
  
  #Try some specific cases:
  #First, try a solution for mostly product = x[0] = [AB] = Ct/2
  cAB = ct/2				# is it your first guess?
  cA=math.sqrt(cAB/KAB)
  cB=cA
  cAp=cA*KA
  cBp=cB*KB
  cA2=cA*cA*KA2
  cB2=cB*cB*KB2
  
  root=fsolve(func,[cAB,cA,cB,cA2,cB2,cAp,cBp])
  print("Trying most product as AB:")
  print(root)
  cdev = ct-(root[0]*2+root[1]+root[2]+2*root[3]+2*root[4]+root[5]+root[6])
  print("concentration deviation = "+str(cdev)+"\n")
  
  
  
  
  #Try the case where unimolecular folding dominates:
  cAp=ct/2
  cBp=ct/2
  cA=cAp/KA
  cB=cBp/KB
  cAB = KAB*cAp*cBp
  cA2=cA*cA*KA2
  cB2=cB*cB*KB2
  
  root=fsolve(func,[cAB,cA,cB,cA2,cB2,cAp,cBp])
  print("Trying most product as unimolecular:")
  print(root)
  cdev = ct-(root[0]*2+root[1]+root[2]+2*root[3]+2*root[4]+root[5]+root[6])
  print("concentration deviation = "+str(cdev)+"\n")
  
  
  #Try the case where homodimer folding dominates:
  cA2=ct/4
  cB2=ct/4
  cA=math.sqrt(cA2/KA2)
  cB=math.sqrt(cB2/KB2)
  cAB = cA*cB*KAB
  cAp=cA*KA
  cBp=cB*KB
  
  
  root=fsolve(func,[cAB,cA,cB,cA2,cB2,cAp,cBp])
  print("Trying most product as homodimer:")
  print(root)
  cdev = ct-(root[0]*2+root[1]+root[2]+2*root[3]+2*root[4]+root[5]+root[6])
  print("concentration deviation = "+str(cdev)+"\n")
  ```