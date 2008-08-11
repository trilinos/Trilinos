import array
from math import *
from sys import *
import Epetra
import time
from Epetra import Comm
from Epetra import SerialComm
from Epetra import Version
from Epetra import Map
from Epetra import Vector
from Epetra import CrsMatrix
from Epetra import CompObject

def cg_method(A, x, b, niters, tol, verbose):
	
	veryVerbose=0
	x.PutScalar(0.0)
	r=Vector(b)
	p=Vector(r)
	w=Vector(A.RowMap())
	rhokm1=r.Dot(r)
	
	counter = A.GetFlopCounter()# DON'T THINK THIS WORKS YET--changed from original once
	
	if(counter !=0):
		p.SetFlopCounter(A)
		r.SetFlopCounter(A)
		w.SetFlopCounter(A)
	
	ierr=1
	betak=0.0
	alphak=0.0
	normr=0.0
	
	
	for iter in range(0,niters):
		if(iter>0):
			betak=rhokm1/rhokm2
			p.Update(1.0,r,betak)
		A.Multiply(false,p,w)
	        ptransw=p.Dot(w)#fix this, you need something different here
		alphak=rhokm1/ptransw
		if(veryVerbose==1):
			print "rhokm1 = " + rhokm1
			print "ptransw = " + ptransw
			print "Alphak = " + alphak
			print "p = " + p 
			print "w = " + w
		x.Update(alphak, p, 1.0)
		r.Update(-alphak, w,1.0)
		rhok=r.Dot(r)
		rhokm2=rhokm1
		rhokm1=rhok
		if(veryVerbose==1):
			print "rhokm2 = " +rhokm2
			print "rhokm1 = " + rhokm1
			print "rhok = "+ rhok
			print "x = " + x
			print "r = " + r
		normr=sqrt(rhok)
		if(iter % 100 == 0 or iter+1==niters or normr < tolerance):
			if(verbose):
				print "Iter = "+iter
				print "Norm of Residual r = b - Ax is " + normr
		if(normr < tolerance):
			ierr=0
			break
	return ierr
def main():
	Comm=Epetra.PyComm()
	MyPID=Comm.MyPID()
	NumProc=Comm.NumProc()
	verbose=(MyPID==0)
	if(verbose):
		print Version()
	print Comm
	
	if(len(argv)!=1): #this needs to be fixed you were being lazy when you pieced this together
		if(verbose): 
			print "Usage: "+argv[0]+" number_of_equations"
		exit(1)
	NumGlobalElements = int(5)
	if (NumGlobalElements < NumProc):
		if (verbose):
			print "numGlobalBlocks = " + NumGlobalElements
			print "cannot be < number of processors = " + NumProc
		exit(1)
	
	PMap=Map(NumGlobalElements, 0, Comm)
	
	NumMyElements = PMap.NumMyElements()
	
	MyGlobalElements = array.array('i', [0]*NumMyElements)
	PMap.MyGlobalElements()
	
	NumNz = array.array('i', [0]*NumMyElements)
	
	for i in range(0,NumMyElements):
		if (MyGlobalElements[i]==0 or MyGlobalElements[i] == NumGlobalElements-1):
			NumNz[i]=2
		else:
			NumNz[i]=3
			
	x=Vector(PMap)
  	b=Vector(PMap)
	b.PutScalar(1.0) 
  	A=CrsMatrix(Epetra.Copy, PMap, 3)#3 should be NumNz
	Values=array.array("d", range(NumMyElements))
	Values[0] = -1.0 
	Values[1] = -1.0
	Indices=array.array("d",range(NumMyElements))
	two=array.array("d")
	two.append(2.0)

	MyGlobalElements2=array.array('i', [0]*NumMyElements)
	
	NumEntries=0

	for i in range(0,NumMyElements):
		if (MyGlobalElements[i]==0):
			Indices[0] = 1
			NumEntries = 1
		elif (MyGlobalElements[i] == NumGlobalElements-1):
			Indices[0] = NumGlobalElements-2
			NumEntries = 1
		else:
			Indices[0] = MyGlobalElements[i]-1
			Indices[1] = MyGlobalElements[i]+1
			NumEntries = 2
		A.InsertGlobalValues(MyGlobalElements[i], Values, Indices)
		A.InsertGlobalValues(MyGlobalElements[i], two, MyGlobalElements2)
		
	ierr = A.FillComplete()
  	assert(ierr==0)
	
	lambd= 0.0
  	niters = NumGlobalElements*10
  	tolerance = .0001
	
	counter =CompObject()
	A.SetFlopCounter(counter)
	start=time.time()#my added lines
	#timer=Epetra_Time(Comm)
	ierr += cg_method(A, x, b, niters, tolerance, verbose)
	#elapsed_time = timer.ElapsedTime()
	end=time.time()#my added lines
	elapsed_time=end-start#my added lines
	total_flops =counter.Flops()
	MFLOPs = total_flops/elapsed_time/1000000.0

  	if (verbose): 
		print "\n\nTotal MFLOPs for first solve = "+MFLOPs
	
	if (verbose): 
		print "\nIncreasing magnitude of first diagonal term, solving again\n\n"
	
	if (A.MyGlobalRow(0)):
		numvals = A.NumGlobalEntries(0)
    		Rowvals =array.array("d",range(numvals))
    		Rowinds =array.array("i" , range(numvals))
    		A.ExtractGlobalRowCopy(0, numvals, numvals, Rowvals, Rowinds)
    		for i in range(0,numvals):
			if (Rowinds[i] == 0): 
				Rowvals[i] *= 10.0

    		A.ReplaceGlobalValues(0, numvals, Rowvals, Rowinds)
	
	lambd = 0.0
  	timer.ResetStartTime()
  	counter.ResetFlops()
  	ierr += cg_method(A, x, b, niters, tolerance, verbose)
  	elapsed_time = timer.ElapsedTime()
  	total_flops = counter.Flops()
  	MFLOPs = total_flops/elapsed_time/1000000.0

  	if (verbose): 
  		print "\n\nTotal MFLOPs for second solve = "+MFLOPs 
		
			
  	MPI_Finalize() 

	return ierr ;
main()
