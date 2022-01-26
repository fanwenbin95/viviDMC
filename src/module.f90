module constants
    implicit none
    
    ! system precision
    integer, parameter :: dp = 8 ! double precision
    
    ! math constants
    real(dp), parameter :: pi = dacos(-1d0)

    ! internal constants
    real(dp), parameter :: fd = 1d-6 ! finite-displacement / delta, unit is Bohr
    real(dp), parameter :: double_fd = 2d0*fd ! $2\delta$
    real(dp), parameter :: quadru_fd_square = 4d0*fd*fd ! $4\delta^2$

contains
    
end module constants


module para
    use constants
    implicit none
    integer :: Natom, Nframe
    integer :: Nmode

    character(4), allocatable :: eleName(:)
    real(dp), allocatable :: eleMass(:)
    real(dp), allocatable :: massMat(:, :)

    ! initial geometries, 3*Natom*Nframe
    real(dp), allocatable :: q0(:, :, :)
    real(dp), allocatable :: grad(:, :)
    real(dp), allocatable :: hess(:, :)

    integer :: Nbond ! to be discarded ! 2022-01-20 23:39:51 Wenbin FAN @FDU
    
    ! DMC parameters
    integer :: N0
    integer :: maxNwalker
    real(dp) :: dt ! time step
    real(dp) :: alpha ! === 0.5 / dt
    real(dp) :: sigma ! === sqrt(dt/M)
    integer :: Nstep ! number of total steps
    integer :: Neq ! number of equilibrium steps

    ! DMC main
    ! Cartesian coordinates of all walkers : 3 * Natom * maxNwalker
    integer, allocatable :: Nwalker(:)
    ! label the ending number of each group of walkers
    integer, allocatable :: batch(:)
    real(dp), allocatable :: q(:, :, :)
    real(dp), allocatable :: Eref(:) ! shifted energy
    real(dp), allocatable :: w(:) ! weights for walkers

    ! sample method
    integer :: sample_type
    ! for vibrational sampling
    real(dp) :: vib_scale
    integer :: Ngrid
    real(dp), allocatable :: vibrc(:) ! 0:Ngrid
    real(dp), allocatable :: vibpot(:, :) ! Ngrid*Nmode
    real(dp), allocatable :: kmat(:, :), pmat(:, :) ! Ngrid*Ngrid
    real(dp), allocatable :: vibwf(:, :), vibcdf(:, :) ! Ngrid*Nmode
    real(dp), allocatable :: vibrand(:, :) ! Nwalker(0)*Nmode

contains
    
end module para


module units
    use constants
    implicit none
    
    ! physical constants
    real(dp), parameter :: clight = 299792458d0

    ! atomic converter
    real(dp), parameter :: au2ang = 0.529177210903d0 ! Angstrom
    real(dp), parameter :: au2kjm = 2625.49963947983d0 ! kJ per mole
    real(dp), parameter :: au2j = 4.3597447222071d-18 ! Joule
    real(dp), parameter :: au2kcm = 627.509474063056d0 ! kcal per mole
    real(dp), parameter :: au2ev = 27.211386245988d0 ! electon Volt
    real(dp), parameter :: au2wn = 219474.6313632d0 ! wavenumber (cm^-1)

    ! mass converter
    real(dp), parameter :: amu2kg = 1.6605390666d-27 ! amu (Dalton) to kg
    real(dp), parameter :: amu2me = 1822.88848621731 ! to electron mass

    ! length converter
    real(dp), parameter :: ang2m = 1d-10 ! Angstrom to meter

contains
    
end module units


module elements
    use constants
    implicit none
    
    integer, parameter :: supNum = 118 ! number of supported elements
    character(2) :: allName(supNum)=(/"H ", "He", & !1~2
    "Li","Be","B ","C ","N ","O ","F ","Ne", & !3~10
    "Na","Mg","Al","Si","P ","S ","Cl","Ar", & !11~18
    "K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", & !19~36
    "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe", & !37~54
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", & !55~71
    "Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", & !72~86
    "Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr", & !87~103
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"/) !104~118
    character(2) :: allNameUp(supNum)=(/"H ", "HE", & !1~2
    "LI","BE","B ","C ","N ","O ","F ","NE", & !3~10
    "NA","MG","AL","SI","P ","S ","CL","AR", & !11~18
    "K ","CA","SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR", & !19~36
    "RB","SR","Y ","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","TE","I ","XE", & !37~54
    "CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU", & !55~71
    "HF","TA","W ","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN", & !72~86
    "FR","RA","AC","TH","PA","U ","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR", & !87~103
    "RF","DB","SG","BH","HS","MT","DS","RG","CN","NH","FL","MC","LV","TS","OG"/) !104~118
    real(dp) :: allMass(supNum)=(/&
      1.0079407541D0,  4.0026019321D0,  6.9400366029D0,  9.0121830650D0, 10.8110280464D0,&
     12.0107358967D0, 14.0067032114D0, 15.9994049243D0, 18.9984031627D0, 20.1800463805D0,&
     22.9897692820D0, 24.3050516198D0, 26.9815385300D0, 28.0854987057D0, 30.9737619984D0,&
     32.0647874061D0, 35.4529375826D0, 39.9477985636D0, 39.0983009101D0, 40.0780225110D0,&
     44.9559082800D0, 47.8667449627D0, 50.9414650374D0, 51.9961317554D0, 54.9380439100D0,&
     55.8451444339D0, 58.9331942900D0, 58.6933471099D0, 63.5460399458D0, 65.3777825295D0,&
     69.7230660726D0, 72.6275501647D0, 74.9215945700D0, 78.9593885570D0, 79.9035277805D0,&
     83.7979999953D0, 85.4676635956D0, 87.6166444696D0, 88.9058403000D0, 91.2236415971D0,&
     92.9063730000D0, 95.9597885412D0, 97.9072124000D0,101.0649401392D0,102.9054980000D0,&
    106.4153275073D0,107.8681496346D0,112.4115578183D0,114.8180866294D0,118.7101125930D0,&
    121.7597836735D0,127.6031264847D0,126.9044719000D0,131.2927614478D0,132.9054519610D0,&
    137.3268916286D0,138.9054688737D0,140.1157307379D0,140.9076576000D0,144.2415960318D0,&
    144.9127559000D0,150.3663557119D0,151.9643781264D0,157.2521306469D0,158.9253547000D0,&
    162.4994728194D0,164.9303288000D0,167.2590826497D0,168.9342179000D0,173.0541501663D0,&
    174.9668149579D0,178.4849787234D0,180.9478756362D0,183.8417775505D0,186.2067045456D0,&
    190.2248596282D0,192.2160516521D0,195.0844568649D0,196.9665687900D0,200.5991670346D0,&
    204.3834128394D0,207.2169080630D0,208.9803991000D0,208.9824308000D0,209.9871479000D0,&
    222.0175782000D0,223.0197360000D0,226.0254103000D0,227.0277523000D0,232.0380558000D0,&
    231.0358842000D0,238.0289104617D0,237.0481736000D0,244.0642053000D0,243.0613813000D0,&
    247.0703541000D0,247.0703073000D0,251.0795886000D0,252.0829800000D0,257.0951061000D0,&
    258.0984315000D0,259.1010300000D0,266.1198300000D0,267.1217900000D0,268.1256700000D0,&
    271.1339300000D0,270.1333600000D0,269.1337500000D0,278.1563100000D0,281.1645100000D0,&
    282.1691200000D0,285.1771200000D0,286.1822100000D0,289.1904200000D0,289.1936300000D0,&
    293.2044900000D0,294.2104600000D0,294.2139200000D0/)

contains
    
end module elements

