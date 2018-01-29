% Blake Leonard	2009

% University of Missouri
% Computational Physics

% balle - Program to compute the trajectory of a baseball using the Euler-Cromer Method

clear;

disp('  balle - Program to compute the trajectory of a baseball');

disp('          using the Euler-Cromer method.');

disp('');


%* Set intitial position and velocity of baseball

y1 = input ('Enter initial height (meters): ');

r1 = [0, y1];       % initial vector position

speed = input ('Enter inital speed (m/s): ');

theta = input ('Enter inital angle (degrees): ');

v1 = [speed*cos(theta*pi/180), speed*sin(theta*pi/180)];   % inital vector velocity

r = r1; v = v1;    % Set initial position and velocity


%* Set physical parameters

Cd = 0.35;                        % Drag Coefficient

area = 4.3e-3;                    % Cross-Sectional Area

grav = 9.81;                      % Gravitational Accel

mass = 0.145;                      % Mass of Projectile

airFlag = input('Air resistance? (Yes:1, No:0): ');

if ( airFlag == 0 )

	rho = 0;                  % No Air Resistance

else

	rho = 1.2;                % Density of Air

end

air_const = -0.5*Cd*rho*area/mass;   % Air Resistance constant


%* Loop until ball hits ground or max steps completed

tau = input('Enter timestep, tau(sec): ');

maxstep = 1000; 

prevx = [0 0 0]; prevy = [0 0 0]; prevt = [0 0 0];    


for jstep = 1:2                     								% 1st step: normal    2nd step: interpolated


	for istep = 1:maxstep


		%* Record position for plotting

		xplot(istep) = r(1);

		yplot(istep) = r(2);

		t = (istep-1)*tau;                 % Current Time

		xNoAir(istep) = r1(1) + v1(1)*t;

		yNoAir(istep) = r1(2) + v1(2)*t - 0.5*grav*t^2;


		%* Calculate the acceleration of the ball

		accel = air_const*norm(v)*v;       % Air Resistance

		accel(2) = accel(2)-grav;          % Gravity
		

		% Interpolate value for r

		if (jstep == 2 && istep > 3)
	
			prevx = [ xplot(istep-3), xplot(istep-2), xplot(istep-1) ];

			prevy = [ yplot(istep-3), yplot(istep-2), yplot(istep-1) ];

			prevt = [ (istep-4)*tau, (istep-3)*tau, (istep-2)*tau ];	

			r(1) = (t-prevt(2))*(t-prevt(3))/((prevt(1)-prevt(2))*(prevt(1)-prevt(3)))*prevx(1) + (t-prevt(1))*(t-prevt(3))/((prevt(2)-prevt(1))*(prevt(2)-prevt(3)))*prevx(2) + (t-prevt(1))*(t-prevt(2))/((prevt(3)-prevt(1))*(prevt(3)-prevt(2)))*prevx(3);

			r(2) = (t-prevt(2))*(t-prevt(3))/((prevt(1)-prevt(2))*(prevt(1)-prevt(3)))*prevy(1) + (t-prevt(1))*(t-prevt(3))/((prevt(2)-prevt(1))*(prevt(2)-prevt(3)))*prevy(2) + (t-prevt(1))*(t-prevt(2))/((prevt(3)-prevt(1))*(prevt(3)-prevt(2)))*prevy(3);

		end
			
	
		%* Calculate the new position and velocity using Euler-Cromer method

		v = v + tau*accel;

		r = r + tau*v;


		%* If ball reaches ground (y<0), break out of the loop

		if ( r(2) < 0 )
	
			xplot(istep+1) = r(1);      % Record last values computed

			yplot(istep+1) = r(2);      

			break;                

		end

	end


	% Compute theoretical No-Air range
	
	for kstep = 1:maxstep						

			t = (kstep-1)*tau;                 % Current Time

			xTheory(kstep) = r1(1) + v1(1)*t;

			yTheory(kstep) = r1(2) + v1(2)*t - 0.5*grav*t^2;

			if ( yTheory(kstep) < 0 )

				theorange = xTheory(kstep);

				break;
			end

	end
	

	%* Print maximum range and time of flight
	
	if ( jstep == 1)

		fprintf('Maximum range is %g meters\n', r(1));

		fprintf('Theoretical No-Air Maximum Range is %g meters\n', theorange);

		origranerror = 100*abs(theorange-r(1))/theorange;

		fprintf('Range Error: %g Percent\n', origranerror );	

		fprintf('Time of flight is %g seconds\n', istep*tau);

		fprintf('Theoretical No-Air Time of flight is %g seconds\n', kstep*tau);

		origtimeerror = 100*abs(kstep*tau-istep*tau)/(kstep*tau);		

		fprintf('Time Error: %g Percent\n', origtimeerror);

		r = r1; v = v1;                    % reset initial position and velocity

	else
		
		fprintf('Corrected Maximum range is %g meters\n', r(1));

		fprintf('Theoretical No-Air Maximum Range is %g meters\n', theorange);

		fprintf('Range Error: %g Percent\n', 100*abs(theorange-r(1))/theorange );	

		fprintf('Range Error Improvement: %g Percent\n', origranerror-100*abs(theorange-r(1))/theorange);

		fprintf('Time of flight is %g seconds\n', istep*tau);

		fprintf('Theoretical No-Air Time of flight is %g seconds\n', kstep*tau);		

		fprintf('Time Error: %g Percent\n', 100*abs(kstep*tau-istep*tau)/(kstep*tau));

		fprintf('Time Error Improvement: %g Percent\n', origtimeerror-100*abs(kstep*tau-istep*tau)/(kstep*tau));
				

	end 


	%* Graph the trajectory of the baseball

	clf; figure(gcf);       % Clear figure window and bring it forward


	% Mark the location of the ground by a straight line

	xground = [0 max(xNoAir)];  yground = [0 0];


	% Plot the computed trajectory and parabolic, no-air curve

	plot(xplot, yplot, '+', xNoAir, yNoAir, '-', xground, yground, '-');
	
	if ( jstep == 1 )

		legend('Euler-Cromer method', 'Theory (No air)', 'Ground' );

	else

		legend('Interpolated Euler-Cromer method', 'Theory (No air)', 'Ground' );

	end
	
	xlabel('Range (m)');  ylabel('Height (m)');

	title('Projectile motion');

	disp('');

	z = input('Hit Enter to Advance');

	disp('');

end