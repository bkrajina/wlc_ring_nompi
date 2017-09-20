#This script evalutes the Alexander polynomial of a closed space curve
#evaluated at t=-1. 

#7/09/2014
#Brad Krajina



from numpy import *
from scipy import *

#load position vector from text file


def alexander(file,folder=''):

    r=loadtxt(file)
    N=len(r[:,1])
    rp=array(zeros((N,3)))  #position vector of projection
    n=array([[0],[0],[1]]) #normal vector of projection plane


    #calculate projection vector

    for i in range(N):
        rp[i,:]=r[i,:]-dot(r[i,:],n)*transpose(n)

    Ncross=0 #Total number of crossings


    Cross=zeros([1,6])

    #Find all crossings

    for i in range(N):
        Ndegen=0 #number of crosses for segment i
    
        for j in range(N):
            #only perform calculations for non-adjacent segments
            if i!=j and i!=j+1 and i!=j-1:
                ip1=i+1
                jp1=j+1
                if i==N-1:
                    ip1=0
                if j==N-1:
                    jp1=0
                    
                if ip1==j or jp1==i:
                    continue
                dri=r[ip1,:]-r[i,:]
                drj=r[jp1,:]-r[j,:]
                #calculate length of segments in the projection and the tangents
                smax=sqrt(sum((rp[ip1,:]-rp[i,:])**2))
                tmax=sqrt(sum((rp[jp1,:]-rp[j,:])**2))
                ui=(rp[ip1,:]-rp[i,:])/smax
                uj=(rp[jp1,:]-rp[j,:])/tmax
                
                udot=dot(ui,uj)
                #if segments are parallel, continue to next segment
                if udot==1. or udot==-1. :
                    continue
                #compute the point of intersection
                tint=(rp[j,1]-rp[i,1]-(ui[1]/ui[0])*(rp[j,0]-rp[i,0]))/((ui[1]*uj[0]/ui[0])-uj[1]);
                sint=(rp[j,0]-rp[i,0]+uj[0]*tint)/ui[0];
                
                #if intersection point is within the length of the line segments
                #count as an intersection
            
                if sint<smax and sint>=0. and tint<tmax and tint>=0.:
                    #determine if this is an over-crossing or under-crossing
                
                    #compute lengths of true segments (not projected)
                    srmax=sqrt(sum((dri**2)))
                    trmax=sqrt(sum((drj**2)))
                    uri=dri/srmax
                    urj=drj/trmax
                    #calculate angle between real vector and projection
                    thetai=arctan(dri[2]/smax)
                    thetaj=arctan(drj[2]/tmax)
                    #calculate the true length along segment of intersection
                    srint=sint/cos(thetai)
                    trint=tint/cos(thetaj)
                    
                    #if this is an undercrossing, storing it in cross
                    if r[i,2]+uri[2]*srint<r[j,2]+urj[2]*trint:
                        Ncross=Ncross+1
                        Cross.resize([Ncross,6])
                        Ndegen=Ndegen+1;
                        Cross[Ncross-1,0]=i;
                        Cross[Ncross-1,1]=j;
                        Cross[Ncross-1,2]=sint;
                        Cross[Ncross-1,3]=tint;
                        Cross[Ncross-1,4]=Ndegen;
                        Cross[Ncross-1,5]=-uj[1]*ui[0]+ui[1]*uj[0]; #cross(ui,uj), + for RH, - for LH

        #after summing over j, set Ndegen for instances
        Cross[Ncross-Ndegen:Ncross,4]=Ndegen



    #sort crosses on same segment according to order of occurence (w.r.t. sint)

    index=0

    while index<Ncross-1:
        Ndegen=Cross[index,4]
        Cross_Slice=copy(Cross[index:index+Ndegen,:])
        Sort_Indices=argsort(Cross_Slice[:,2])
        Cross[index:index+Ndegen,:]=Cross_Slice[Sort_Indices,:]
        index=index+Ndegen
   

    #Construct vector of over-pass indices according to cross indexing described by Vologodskii

    over_ind=zeros([Ncross,1])
    over_ind.dtype=int
    #Loop over all crossings
       
    for i in range(Ncross):
        j=Cross[i,1] #j index of over-pass
    
    
        #special cases. j comes before the first crossing or after the last crossing
        if j<Cross[0,0]:
            over_ind[i]=0
            continue
        elif j>Cross[Ncross-1,0]:
            over_ind[i]=Ncross
            continue
       
        #sum over all crossings to determine which under-crossings this lies between
        for k in range(Ncross-1):
            #if j lies between cross k and cross k+1, it is segment k+1
            if j>Cross[k,0] and j<Cross[k+1,0]:
                over_ind[i]=k+1
                break
                #if j=i for cross k, need to determine where over-pass is relative to under-passings
            elif j==Cross[k,0]:
                t=Cross[i,3]
                Ndegen=Cross[k,4]
                #special cases: t is before the first under-pass or the last under-pass of segment j
                if t<=Cross[k,2]:
                    over_ind[i]=k
                    break
                elif t>Cross[k+Ndegen-1,2]:
                    over_ind[i]=k+Ndegen
                    break
                    #other-wise, determine which undercrossings t lies between
                else: 
                    ind=1
                    #loop over all degenerate undercrossings
                    while ind<Ndegen:
                        #if t lies between the s of k+ind-1 undercrossing and the next,
                        #then this overpass has a new index of k+ind
                        if t>Cross[k+ind-1,2] and t<=Cross[k+ind,2]:
                            over_ind[i]=k+ind
                            break
                        ind=ind+1    
                    break

    #Calculate the Alexander matrix
    #Evaluated at t=-1


    A=zeros([Ncross,Ncross])

    for k in range(Ncross):
        kp1=k+1
        i=over_ind[k]
        #periodic index conditions
        if k==Ncross-1:
            kp1=0
        if i==Ncross:
            i=0
            #calculate values
        if i==k or i==kp1:
            A[k,k]=-1
            A[k,kp1]=1
        else:
            A[k,k]=1
            A[k,kp1]=1
            A[k,i]=-2

    #calculate determinant of matrix with one row and one column removed

    if A[0:Ncross-1,0:Ncross-1].size==0:
        #if A has one crossing or less, trivial knot
        delta=1
    else:
        delta=linalg.det(A[0:Ncross-1,0:Ncross-1])

    fileind=file.split('r')[-1]    
#    if folder!='':
 #       savetxt(folder+'/cross'+fileind,Cross)
  #      savetxt(folder+'/over_ind_'+fileind,over_ind)
   #     savetxt(folder+'/A_'+fileind,A)
                

    return delta    
        
    
    
                           
                
        
                
            





