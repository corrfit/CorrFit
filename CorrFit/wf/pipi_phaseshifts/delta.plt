clear
orientation portrait
device\colour postscript
!! cursor -2 is for centering, -1 is for left justify, -3 for right just.
set cursor -2

!___________________________________________
set xleadz 1.0
set yleadz 1.0
set xnumsz .6
set ynumsz .6
set xticl -.5
set xtics -.25
set yticl -.5
set ytics -.25
set xlabsz 1.0
set ylabsz 1.0
set lintyp 1.0
set linthk 6.0
set charsz 0.4
!set xlog 10
!set ylog 10
!!ylabel "P(n)"
!!xlabel "n"
!___________________________________________
!These are the dimensions of the borders
xxll=4.5
xxuu=17
yyll=5.0
yyuu=15
!___________________________________________
!delxlabel and delylabel are the offsets of the label.
delxlabel=3
delylabel=3.5
set xlaxis xxll
set xuaxis xxuu
set ylaxis yyll
set yuaxis yyuu
!___________________________________________
color 1

! y-axis label
ytext=(yyll+yyuu)/2.0
xtext=xxll-delylabel
set xloc xtext
set yloc ytext
set txtang 90
set txthit 0.7
text `<delta> (degrees)'

! x-axis label
xtext=(xxuu+xxll)/2.0
ytext=yyll-delxlabel
set xloc xtext
set yloc ytext
set txtang 0
set txthit .7
text `E (MeV/c)'

!autoscale
!usage: scales xmin,xmax,ymin,ymax
scales 0 5000 0 180

! These 8 par.s set the locations of tick marks
set xvmin 0.0
set xvmax 2000.0
set nlxinc 4
set nsxinc 5

set yvmin -180
set yvmax 180
set nlyinc 4
set nsyinc 9

graph\axesonly

vector q p 80
n=0
do I = [1:80:1]
q(I)=I*50
enddo
e=2*sqrt(q*q+139.58*139.58)

a1=0.214/197.323
a2=290.0
a3=625
p=(180.0/3.141592654)*atan(a1*q/(1.0-(q/a2)+q*q/(a3*a3)))
p=p+180.0*(abs(p)-p)/abs(2*p)

color red
set pchar 12
graph\noaxes e p

vector x y 2
x(1)=0
x(2)=2000
set pchar 0
y(1)=90
y(2)=90
graph\noaxes x y
y(1)=45
y(2)=45
graph\noaxes x y

!-----------------------------------------------

hardcopy s delta.eps
