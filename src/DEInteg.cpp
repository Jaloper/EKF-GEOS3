#include "..\include\DEInteg.hpp"
#include "..\include\sign_.hpp"


#define eps 2.2204e-16

Matrix DEInteg(Matrix func(double t, Matrix y), double t,double tout,double relerr,double abserr,int n_eqn,Matrix& y){

	// maxnum = 500;
	double twou  = 2*eps;
	double fouru = 4*eps;

	enum class DE_STATE {
		DE_INIT = 1,      // Restart integration
		DE_DONE = 2,      // Successful step
		DE_BADACC = 3,    // Accuracy requirement could not be achieved
		DE_NUMSTEPS = 4,  // Permitted number of steps exceeded
		DE_STIFF = 5,     // Stiff problem suspected
		DE_INVPARAM = 6   // Invalid input parameters
	};
	
	DE_STATE State_ = DE_STATE::DE_INIT;
	bool PermitTOUT = true;         // Allow integration past tout by default
	double told = 0.0;

	// Powers of two (two(n)=2^n)
	Matrix two(14);  
	two(1)=1.0; two(2)=2.0; two(3)=4.0; two(4)=8.0; two(5)=16.0; two(6)=32.0; two(7)=64.0; two(8)=128.0;
			two(9)=256.0; two(10)=512.0; two(11)=1024.0; two(12)=2048.0; two(13)=4096.0; two(14)=8192.0;

	Matrix gstr(14);
	gstr(1)=1.0; gstr(2)= 0.5; gstr(3)=0.0833; gstr(4)=0.0417; gstr(5)=0.0264; gstr(6)=0.0188; gstr(7)=0.0143; gstr(8)=0.0114;
			gstr(9)=0.00936; gstr(10)=0.00789; gstr(11)=0.00679; gstr(12)=0.00592; gstr(13)=0.00524; gstr(14)=0.00468;

	Matrix yy    = zeros(1,n_eqn);    // Allocate vectors with proper dimension
	Matrix wt    = zeros(1,n_eqn);
	Matrix p     = zeros(1,n_eqn);
	Matrix yp    = zeros(n_eqn,1);
	Matrix phi   = zeros(n_eqn,17);
	Matrix g     = zeros(1,14);
	Matrix sig   = zeros(1,14);
	Matrix rho   = zeros(1,14);
	Matrix w     = zeros(1,13);
	Matrix alpha = zeros(1,13);
	Matrix beta  = zeros(1,13);
	Matrix v     = zeros(1,13);
	Matrix psi_  = zeros(1,13);

	// while(true)

	// Return, if output time equals input time

	if (t==tout)    // No integration
		return y;

	// Test for improper parameters

	double epsilon = max(relerr,abserr);

	if ( ( relerr <  0.0) ||        // Negative relative error bound
		 ( abserr <  0.0         ) ||        // Negative absolute error bound
		 ( epsilon    <= 0.0         ) ||     // Both error bounds are non-positive
		 ( State_  >  DE_STATE::DE_INVPARAM ) ||  // Invalid status flag
		 ( (State_ != DE_STATE::DE_INIT) &&    
		 (t != told)           ) ){
			 State_ = DE_STATE::DE_INVPARAM;              // Set error code
			 return y;                                  // Exit
		 }

	// On each call set interval of integration and counter for
	// number of steps. Adjust input error tolerances to define
	// weight vector for subroutine STEP.

	double del    = tout - t;
	double absdel = abs(del);

	double tend   = t + 100.0*del;
	if (!PermitTOUT)	tend = tout;

	int nostep = 0;
	int kle4   = 0;
	bool stiff  = false;
	double releps = relerr/epsilon;
	double abseps = abserr/epsilon;
	double h, x;
	double delsgn;
	bool OldPermit, start;
	if  ( (State_==DE_STATE::DE_INIT) || (!OldPermit) || (delsgn*del<=0.0) ){
		// On start and restart also set the work variables x and yy(*),
		// store the direction of integration and initialize the step size
		 start  = true;
		x      = t;
		yy     = transpose(y);
		delsgn = std::copysign(1.0, del);
		h      = std::copysign( max(fouru*abs(x), abs(tout-x)), tout-x );
	}
	double temp1;
	double hold;
	 double hnew;
	 bool phase1;
	 bool nornd;
	 double absh;
	 int k, kold;
	while (true){   // Start step loop

	  // If already past output point, interpolate solution and return
	  if (abs(x-t) >= absdel){
		  Matrix yout  = zeros(n_eqn,1);
		  Matrix ypout = zeros(n_eqn,1);
		  g(2)= 1.0;
		  rho(2)= 1.0;
		  double hi = tout - x;
		  int ki = kold + 1;
		  
		  // Initialize w[*] for computing g[*]
		  Matrix w=zeros(1,ki);
		  for(int i = 1; i <= ki; i++){
			  temp1 = i;
			  w(i+1) = 1.0/temp1;
			}
		  // Compute g[*]
		  double term = 0.0;
		  for(int j = 2; j <= ki; j++){
			  double psijm1 = psi_(j);
			  double gamma = (hi + term)/psijm1;
			  double eta = hi/psijm1;
			  for(int i = 1; i <= ki+1-j; i++){
				  w(i+1) = gamma*w(i+1) - eta*w(i+2);
			  }
			  g(j+1) = w(2);
			  rho(j+1) = gamma*rho(j);
			  term = psijm1;
		  }
		  
		  // Interpolate for the solution yout and for
		  // the derivative of the solution ypout 
			yout=zeros(n_eqn,1);
		  for(int j = 1; j <= ki; j++){
			  int i = ki+1-j;
			  yout  = yout  + extract_column(phi,i+1)*g(i+1);
			  Matrix ypout = ypout + extract_column(phi,i+1)*rho(i+1);
		  }
		  yout = y + yout*hi;
		  y    = yout;
		  State_    = DE_STATE::DE_DONE; // Set return code
		  t         = tout;             // Set independent variable
		  told      = t;                // Store independent variable
		  OldPermit = PermitTOUT;
		  return y;                       // Normal exit
		}   
		
	  // If cannot go past output point and sufficiently close,
	  // extrapolate and return
	  if ( !PermitTOUT && ( abs(tout-x) < fouru*abs(x) ) ){
		  h = tout - x;
		  yp = func(x,yy);          // Compute derivative yp(x)
		  y = yy + yp*h;                // Extrapolate vector from x to tout
		  State_    = DE_STATE::DE_DONE; // Set return code
		  t         = tout;             // Set independent variable
		  told      = t;                // Store independent variable
		  OldPermit = PermitTOUT;
		  return y;                       // Normal exit
	  }
	  
	  // Test for too much work
	//   if (nostep >= maxnum)
	//       State_ = DE_STATE::DE_NUMSTEPS; // Too many steps
	//       if (stiff) 
	//           State_ = DE_STATE::DE_STIFF;// Stiffness suspected
	//       end
	//       y         = yy;                // Copy last step
	//       t         = x;
	//       told      = t;
	//       OldPermit = true;
	//       return;                        // Weak failure exit
	//   end
	  
	  // Limit step size, set weight vector and take a step
	  h  = sign_(min(abs(h), abs(tend-x)), h);
	  for (int l=1; l<=n_eqn;l++){
		  wt(l) = releps*abs(yy(l)) + abseps;
	  }
		std::cout <<"BBBBBBBBBBBBBBBBBBBBBBBBBBB\n"<<endl;
	//   Step
	//                                                                   
	// Begin block 0                                                     
	//                                                                   
	// Check if step size or error tolerance is too small for machine    
	// precision.  If first step, initialize phi array and estimate a    
	// starting step size. If step size is too small, determine an       
	// acceptable one.                                                   
	//                                                                   
	bool crash;
	if (abs(h) < fouru*abs(x)){
		h = sign_(fouru*abs(x),h);
		crash = true;
		return y;           // Exit 
	}
	
	double p5eps  = 0.5*epsilon;
	crash  = false;
	g(2)   = 1.0;
	g(3)   = 0.5;
	sig(2) = 1.0;
	int ifail = 0;
		std::cout <<"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"<<endl;
	// If error tolerance is too small, increase it to an 
	// acceptable value.                                  

	double round = 0.0;
	y=transpose(y);
	for (int l=1; l<=n_eqn;l++){
		round = round + (y(l)*y(l))/(wt(l)*wt(l));
	}
	y=transpose(y);
	
	round = twou*sqrt(round);
	if (p5eps<round){
		epsilon = 2.0*round*(1.0+fouru);
		crash = true;
		return y;
	}
		std::cout <<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\n"<<endl;
	if (start){
	  // Initialize. Compute appropriate step size for first step. 
	  yp = func(x,y);
	  yp=transpose(yp);
	  double sum = 0.0;
	  for (int l=1; l<=n_eqn;l++){
		  phi(l,2) = yp(l);
		  phi(l,3) = 0.0;
		  sum = sum + (yp(l)*yp(l))/(wt(l)*wt(l));
	  }
	  sum  = sqrt(sum);
	  absh = abs(h);
	  if (epsilon<16.0*sum*h*h){
		  absh=0.25*sqrt(epsilon/sum);
	  }
	  h    = sign_(max(absh, fouru*abs(x)), h);
	  hold = 0.0;
	  hnew = 0.0;
	  k    = 1;
	  kold = 0;
	  start  = false;
	  phase1 = true;
	  nornd  = true;
	  if (p5eps<=100.0*round){
		  nornd = false;
		 for (int l=1; l<=n_eqn;l++){
			  phi(l,16)=0.0;
		 }
	  }
	}
		std::cout <<"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n"<<endl;
	//                                                                   
	// End block 0                                                       
	//                                                                   

	//                                                                   
	// Repeat blocks 1, 2 (and 3) until step is successful               
	//        
		int kp1;
	  int kp2;
	  int km1;
	  int km2;	
	  int knew;
	  int ns;
	  int nsp1,nsp2;
	  int realns;
	  int im1;
	  double temp2,temp3, temp4, temp5, temp6;
	  int reali;
	  int nsm2, it;
	  double erk,erkm1,erkm2;
	while(true){
		
	  //                                                                 
	  // Begin block 1                                                   
	  //                                                                 
	  // Compute coefficients of formulas for this step. Avoid computing 
	  // those quantities not changed when step size is not changed.     
	  //                                                                 
	  
	  kp1 = k+1;
	  kp2 = k+2;
	  km1 = k-1;
	  km2 = k-2;
	  
	  // ns is the number of steps taken with size h, including the 
	  // current one. When k<ns, no coefficients change.      
	  
	  if (h !=hold)
		  ns=0;
	  if (ns<=kold)
		  ns=ns+1;
	  nsp1 = ns+1;
			std::cout <<"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"<<endl;
	  if (k>=ns){
		  // Compute those components of alpha[*],beta[*],psi[*],sig[*] 
		  // which are changed                                          
		  beta(ns+1) = 1.0;
		  realns = ns;
		  alpha(ns+1) = 1.0/realns;
		  temp1 = h*realns;
		  sig(nsp1+1) = 1.0;
		  if (k>=nsp1){
			  for (int i=nsp1;i<=k;i++){
				  im1   = i-1;
				  temp2 = psi_(im1+1);
				  psi_(im1+1) = temp1;
				  beta(i+1)  = beta(im1+1)*psi_(im1+1)/temp2;
				  temp1    = temp2 + h;
				  alpha(i+1) = h/temp1;
				  reali = i;
				  sig(i+2) = reali*alpha(i+1)*sig(i+1);
			  }
		  }
		  psi_(k+1) = temp1;
		  // Compute coefficients g[*]; initialize v[*] and set w[*].
		  if (ns>1){
			  // If order was raised, update diagonal part of v[*]
			  if (k>kold){
				  temp4 = k*kp1;
				  v(k+1) = 1.0/temp4;
				  nsm2 = ns-2;
				  for (int j=1;j<=nsm2;j++){
					  it = k-j;
					  v(it+1) = v(it+1) - alpha(j+2)*v(it+2);
				  }
			  }
			  
			  // Update V[*] and set W[*]
			  int limit1 = kp1 - ns;
			  temp5  = alpha(ns+1);
			  for (int iq=1;iq<limit1;iq++){
				  v(iq+1) = v(iq+1) - temp5*v(iq+2);
				  w(iq+1) = v(iq+1);
			  }
			  g(nsp1+1) = w(2);
		  }
		  else{
			  for (int iq=1;iq<=k;iq++){
				  temp3 = iq*(iq+1);
				  v(iq+1) = 1.0/temp3;
				  w(iq+1) = 1.0/temp3;
			  }
		  }
		std::cout <<"HHHHHHHHHHHHHHHHHHHHHHHH\n"<<endl;
		  // Compute the g[*] in the work vector w[*]
		  nsp2 = ns + 2;
		  if (kp1>=nsp2){
			  for(int i=nsp2;i<=kp1;i++){
				  int limit2 = kp2 - i;
				  temp6  = alpha(i);
				  for (int iq=1;iq<=limit2;iq++){
					  w(iq+1) = w(iq+1) - temp6*w(iq+2);
				  }
				  g(i+1) = w(2);
			  }
		  }
		  
	  } // if K>=NS
	  
	  //
	  // End block 1
	  //
	  
	  //
	  // Begin block 2
	  //
	  // Predict a solution p[*], evaluate derivatives using predicted
	  // solution, estimate local error at order k and errors at orders
	  // k, k-1, k-2 as if constant step size were used.
	  //   
	  
	  // Change phi to phi star
	  if (k>=nsp1){
		  for (int i=nsp1;i<=k;i++){
			  temp1 = beta(i+1);
			  for (int l=1; l<n_eqn;l++){
				  phi(l,i+1) = temp1 * phi(l,i+1);
			  }
		  }
	}
		std::cout <<"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"<<endl;
	  // Predict solution and differences 
	   for (int l=1; l<=n_eqn;l++){
		  phi(l,kp2+1) = phi(l,kp1+1);
		  phi(l,kp1+1) = 0.0;
		  p(l)       = 0.0;
	   }
	   int ip1;
	  for (int j=1;j<=k;j++){
		  it     = kp1 - j;
		  ip1   = it+1;
		  temp2 = g(it+1);
		  for (int l=1;l<=n_eqn;l++){
			  p(l)     = p(l) + temp2*phi(l,it+1);
			  phi(l,it+1) = phi(l,it+1) + phi(l,ip1+1);
		  }
	  }
		std::cout <<"JJJJJJJJJJJJJJJJJJJJJJJJJ\n"<<endl;
	  if (nornd)
		  p = transpose(y) + p*h;
	  else{
		  for (int l=1;l<=n_eqn;l++){
			  double tau = h*p(l) - phi(l,16);
			  p(l) = y(l) + tau;
			  phi(l,17) = (p(l) - y(l)) - tau;
		  }
	  }
		std::cout <<"KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK\n"<<endl;
	  double xold = x;
	  x = x + h;
	  absh = abs(h);
	  yp = func(x,transpose(p));
	  yp=transpose(yp);
	  
	  // Estimate errors at orders k, k-1, k-2 
	  erkm2 = 0.0;
	  erkm1 = 0.0;
	  erk = 0.0;
		std::cout <<"LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL\n"<<endl;
		
	  for (int l=1;l<=n_eqn;l++){
		  temp3 = 1.0/wt(l);
		  temp4 = yp(l) - phi(l,1+1);
		   //std::cout <<"temp3\n"<< temp3<<endl;
		  std::cout <<"temp4\n"<< temp4<<endl;
		  std::cout <<"yp(l)\n"<< yp(l)<<endl;
		  std::cout <<"phi(l,1+1)\n"<< phi(l,1+1)<<endl;
		  if (km2> 0){
			  erkm2 = erkm2 + ((phi(l,km1+1)+temp4)*temp3)
							 *((phi(l,km1+1)+temp4)*temp3);
		  }
		  if (km2>=0){
			  erkm1 = erkm1 + ((phi(l,k+1)+temp4)*temp3)
							 *((phi(l,k+1)+temp4)*temp3);
		  }
		  erk = erk + (temp4*temp3)*(temp4*temp3);
	  }
	  
	  if (km2> 0)
		  erkm2 = absh*sig(km1+1)*gstr(km2+1)*sqrt(erkm2);
	  
	  if (km2>=0)
		  erkm1 = absh*sig(k+1)*gstr(km1+1)*sqrt(erkm1);
	  
		std::cout <<"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\n"<<endl;
	  temp5 = absh*sqrt(erk);
	 double err = temp5*(g(k+1)-g(kp1+1));
	  erk = temp5*sig(kp1+1)*gstr(k+1);
	 knew = k;
	  //std::cout <<"km2\n"<< km2<<endl;
	  //std::cout <<"erk\n"<< erk<<endl;
	  //std::cout <<"erkm1\n"<< erkm1<<endl;
	  //std::cout <<"erkm2\n"<< erkm2<<endl;
	  // Test if order should be lowered 
	  if (km2 >0){
		  if (max(erkm1,erkm2)<=erk){
			  std::cout <<"AAAAAA\n"<<endl;
			  knew=km1;
		  }
	  }
	  if (km2==0){
		  if (erkm1<=0.5*erk){
			  knew=km1;
		  }
	  }
	  
	  //
	  // End block 2
	  //
	  
	  //
	  // If step is successful continue with block 4, otherwise repeat
	  // blocks 1 and 2 after executing block 3
	  //
	  
	  bool success = (err<=epsilon);
		std::cout <<"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"<<endl;
	  if (!success){	  
		//
		// Begin block 3
		//
		
		// The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
		// 3rd consecutive failure, set order to 1. If step fails more
		// than 3 times, consider an optimal step size. Double error
		// tolerance and return if estimated step size is too small
		// for machine precision.
		//
		
		// Restore x, phi[*,*] and psi[*]
		phase1 = false; 
		x = xold;
		for (int i=1;i<=k;i++){
			temp1 = 1.0/beta(i+1);
			ip1 = i+1;
			for (int l=1;l<=n_eqn;l++){
				phi(l,i+1)=temp1*(phi(l,i+1)-phi(l,ip1+1));
			}
		}
		
		if (k>=2){
			for (int i=2;i<=k;i++){
				psi_(i) = psi_(i+1) - h;
			}
		}
		std::cout <<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n"<<endl;
		// On third failure, set order to one. 
		// Thereafter, use optimal step size   
		ifail = ifail+1;
		temp2 = 0.5;
		if (ifail>3){ 
		  if (p5eps < 0.25*erk){
			  temp2 = sqrt(p5eps/erk);
		  }
		}
		if (ifail>=3){
			knew = 1;
		}
		h = temp2*h;
		k = knew;
		if (abs(h)<fouru*abs(x)){
			crash = true;
			h = sign_(fouru*abs(x), h);
			epsilon = epsilon*2.0;
			return y;                 // Exit 
		}
		
		//
		// End block 3, return to start of block 1
		//
		
	  }  // end if(success)
		std::cout <<"PPPPPPPPPPPPPPPPPPPPPPPPPPPPPP\n"<<endl;
		std::cout <<y<<endl;	  
	  if (success)
		  break;
	  
	  
	}
	
	//
	// Begin block 4
	//
	// The step is successful. Correct the predicted solution, evaluate
	// the derivatives using the corrected solution and update the
	// differences. Determine best order and step size for next step.
	//

	kold = k;
	hold = h;
	y=transpose(y);
	// Correct and evaluate
	std::cout << "Valor de h: " << h << std::endl;
std::cout << "Valor de g(kp1+1): " << g(kp1+1) << std::endl;
std::cout << "kp1 " << kp1 << std::endl;
std::cout << "Valores de yp:\n " << yp << std::endl;
std::cout << "Valores de phi(col 2):\n " << extract_column(phi, 2) << std::endl;
	
	temp1 = h*g(kp1+1);
	if (nornd){
		for (int l=1;l<=n_eqn;l++){
			y(l) = p(l) + temp1*(yp(l) - phi(l,2));
		}
	}
	else{
		for (int l=1;l<=n_eqn;l++){
			double rho_aux = temp1*(yp(l) - phi(l,2)) - phi(l,17);
			y(l) =  rho_aux+ p(l);
			phi(l,16) = (y(l) - p(l)) - rho_aux;
		}
	}
		std::cout <<"QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQqq\n"<<endl;
		std::cout <<transpose(y)<<endl;	
	yp = func(x,transpose(y));
	yp=transpose(yp);
	y=transpose(y);
	
	// Update differences for next step 
	for (int l=1;l<=n_eqn;l++){
		phi(l,kp1+1) = yp(l) - phi(l,2);
		phi(l,kp2+1) = phi(l,kp1+1) - phi(l,kp2+1);
	}
	
	for (int i=1;i<=k;i++){
		for (int l=1;l<=n_eqn;l++){
			phi(l,i+1) = phi(l,i+1) + phi(l,kp1+1);
		}
	}
	
	// Estimate error at order k+1 unless               
	// - in first phase when always raise order,        
	// - already decided to lower order,                
	// - step size not constant so estimate unreliable  
	double erkp1 = 0.0;
	//std::cout << "knew " << knew << std::endl;
	//std::cout << "km1 " << km1 << std::endl;
	//std::cout << "phase1 " << phase1 << std::endl;
	//std::cout << "erk " << erk << std::endl;
	//std::cout << "k " << k << std::endl;
	if ( (knew==km1) || (k==12) )
		phase1 = false;

	if (phase1){
		k = kp1;
		erk = erkp1;
	}
	else{
		if (knew==km1){
			//std::cout <<"knew==km1\n"<<endl;
			// lower order 
			k = km1;
			erk = erkm1;
		}
		else{
			if (kp1<=ns){
				//std::cout <<"kp1<=ns\n"<<endl;
				for (int l=1;l<=n_eqn;l++){
					erkp1 = erkp1 + (phi(l,kp2+1)/wt(l))*(phi(l,kp2+1)/wt(l));
				}
				erkp1 = absh*gstr(kp1+1)*sqrt(erkp1);
				// Using estimated error at order k+1, determine 
				// appropriate order for next step               
				if (k>1){
					//std::cout <<"k>1"<<endl;
					if ( erkm1<=min(erk,erkp1)){
						// lower order
						k=km1; erk=erkm1;
					}else{
						if ( (erkp1<erk) && (k!=12) ){
							// raise order 
							k=kp1;
							erk=erkp1;
						}
					}
				}else{ if (erkp1<0.5*erk){
					//std::cout <<"erkp1<0.5*erk\n"<<endl;
					// raise order
					// Here erkp1 < erk < max(erkm1,ermk2) else    
					// order would have been lowered in block 2.   
					// Thus order is to be raised                  
					k = kp1;
					erk = erkp1;
				}}
			} // end if kp1<=ns
		} // end if knew!=km1
	} // end if !phase1
	//std::cout << "k " << k << std::endl;
	//std::cout << "erk " << erk << std::endl;
		std::cout <<"RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR\n"<<endl;
	// With new order determine appropriate step size for next step
	if ( phase1 || (p5eps>=erk*two(k+2)) ){
		hnew = 2.0*h;
	}else{
		if (p5eps<erk){
			temp2 = k+1;
			double r = p5eps / pow(erk, 1.0/temp2);
			hnew = absh*max(0.5, min(0.9,r));
			hnew = sign_(max(hnew, fouru*abs(x)), h);
		}else{
			hnew = h;
		}
	}
	h = hnew;

	//
	// End block 4
	//
	
	  // Test for too small tolerances
	  if (crash){
		  State_    = DE_STATE::DE_BADACC;
		  relerr    = epsilon*releps;       // Modify relative and absolute
		  abserr    = epsilon*abseps;       // accuracy requirements
		  y         = yy;                   // Copy last step
		  t         = x;
		  told      = t;
		  OldPermit = true;
		  return y;                       // Weak failure exit
	  }
	  
	  nostep = nostep+1;  // Count total number of steps
	 
		std::cout <<"SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs\n"<<endl;
		std::cout <<"y\n"<<y<<endl;	
		std::cout <<"nostep\n"<<nostep<<endl;	
																											if(nostep==15) return y;		
	  // Count number of consecutive steps taken with the order of
	  // the method being less or equal to four and test for stiffness
	  kle4 = kle4+1;
	  if (kold>  4){
		  kle4 = 0;
	  }
	  if (kle4>=50){
		  stiff = true;
	  }
	  
	} // End step loop
	  
	//   if ( State_==DE_STATE::DE_INVPARAM )
	//       error ('invalid parameters in DEInteg');
	//       exit; 
	//   end
	//   if ( State_==DE_STATE::DE_BADACC )
	//       warning ('on','Accuracy requirement not achieved in DEInteg');
	//   end
	//   if ( State_==DE_STATE::DE_STIFF )
	//       warning ('on','Stiff problem suspected in DEInteg');
	//   end
	//   if ( State_ >= DE_STATE::DE_DONE )
	//       break;
	//   end
	//   
	// end
return y;
}