%
% Version 21.07.2015  - fitted DOS between symmetrie points
% Version 26.07.2015  - Correct Slater-Koster d-d
% Version 24.08.2015  - Approximated Total Energy
% Version 01.11.2015  - covalent radii according Cordero
% Version 14.11.2015  - metal radii for metals
% Version 19.11.2015  - Optical Properties
% Version 28.12.2015  - Hexagonale structures
% Version 31.05.2017  - approx. f-electrons in DOS
% Version 25.11.2018  - Plot Structure
%
clear all
alat = 0
nkmax = 6
nkbfit = 50

quiet silent
maxcompund = 3
maxelemente = 89
atomlabel =          [  'H '; 'He']
atomlabel = [atomlabel; 'Li'; 'Be'; 'B '; 'C '; 'N '; 'O '; 'F '; 'Ne']
atomlabel = [atomlabel; 'Na'; 'Mg'; 'Al'; 'Si'; 'P '; 'S '; 'Cl'; 'Ar']
atomlabel = [atomlabel; 'K '; 'Ca']
atomlabel = [atomlabel;             'Sc'; 'Ti'; 'V '; 'Cr'; 'Mn'; 'Fe'; 'Co'; 'Ni'; 'Cu'; 'Zn']
atomlabel = [atomlabel;             'Ga'; 'Ge'; 'As'; 'Se'; 'Br'; 'Kr']
atomlabel = [atomlabel; 'Rb'; 'Sr']
atomlabel = [atomlabel;             'Y '; 'Zr'; 'Nb'; 'Mo'; 'Tc'; 'Ru'; 'Rh'; 'Pd'; 'Ag'; 'Cd']
atomlabel = [atomlabel;             'In'; 'Sn'; 'Sb'; 'Te'; 'I '; 'Xe'];
atomlabel = [atomlabel; 'Cs'; 'Ba']
atomlabel = [atomlabel;             'La'; 'Ce'; 'Pr'; 'Nd'; 'Pm'; 'Sm'; 'Eu'; 'Gd'; 'Tb'; 'Dy'; 'Ho'; 'Er'; 'Tm'; 'Yb']
atomlabel = [atomlabel;             'Lu'; 'Hf'; 'Ta'; 'W '; 'Re'; 'Os'; 'Ir'; 'Pt'; 'Au'; 'Hg']
atomlabel = [atomlabel;             'Tl'; 'Pb'; 'Bi'; 'Po'; 'At'; 'Rn'];
atomlabel = [atomlabel; 'LK'; 'KT'; 'AN']
quiet normal

disp 'Klick <Enter> to start' 
dummy = input('Start')
% disp 'If the lattice constant is unknown tpye 0. The lattice constant ist calculated from the listed atomic radii' 
% alat = input('Lattice constant: ')
disp 'Possible Structuretypes are:'
disp 'Ah   - simple cubic phase (alpha-Po)'
disp 'A1   - cubic face centered dense packing (Ca)'
disp 'A2   - cubic body centered packing (Na)'
disp 'A3   - hexagonal closed packing (Mg)'
disp 'A4   - Tetrahedral Structure (Diamond)'
disp 'B1   - face centerd cubic rocksalt structure (NaCl)'
disp 'B2   - simple cubic body center structure (CsCl)'
disp 'B3   - Sphalerit type mixed tetrahedral cubic structure (ZnS)'
disp 'B_Z  - Zintl-Phase cubic structure (NaTl)'
disp 'C1   - tetrahedral cubic face center structure (CaF2)'
disp 'C3   - Cuprit cubic face center structure (Cu2O)'
disp 'C9   - beta-Crystobalit structure (pseudo) fcc (SiO2)'
disp 'C15  - Laves phase cubic structure (MgCu2)'
disp 'E2_1 - Perowskit structure (CaTiO3)'
disp 'E2_X - pseudo AlPO4- cubic face center structure (ABC3)'
disp 'L1_2 - cubic face center structure (Cu3Au)'
disp 'Q_B3 - Quartary tetrahedral compound (AB)(CD) (CuGaTe2)'
disp 'T_C15 - Quasibinäre Laves Phase (A(BC)2 (TiFeAl)'
quiet silent
inStructure = input('Select a Structuretype from the list above: ','s')

if (isequal(inStructure(1),'A'))
    maxcompund=1
end
if (isequal(inStructure(1),'B'))
    maxcompund=2
end
if (isequal(inStructure(1),'C'))
    maxcompund=2
end
if (isequal(inStructure(1),'E'))
    maxcompund=3
end
if (isequal(inStructure(1),'H'))
    maxcompund=3
end
if (isequal(inStructure(1),'L'))
    maxcompund=2
end
if (isequal(inStructure(1),'T'))
    maxcompund=3
end
if (isequal(inStructure(1),'Q'))
    maxcompund=4
end

disp 'Now enter the different compunds step by step:'
inu = 0
for (ic = 1:maxcompund )
    Compound = input('Compound: ','s')
    lc = length(Compound)
    if (lc > 0)   
        atomtyp(ic,1:lc) = Compound
        if (isequal(lc,1))
            atomtyp(ic,2) = ' '
        end
        if (lc < 3)
            for (ie = 1:maxelemente)
                if (isequal(atomtyp(ic,1:2),atomlabel(ie,1:2)))
                    inu = inu + 1
                    atomselect(inu) = ie
                end
            end
        end
    end
end



if isequal(inStructure,'Ah')
    cryst = 'cub'
    Nnu = 1;
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
end

if isequal(inStructure,'Ah_d')
    cryst = 'cub'
    Nnu = 8
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.000 ; 0.000 ; 0.500 ] 
    atomchg(3) = atomselect(1); atompos(1:3,3) = [ 0.000 ; 0.500 ; 0.000 ] 
    atomchg(4) = atomselect(1); atompos(1:3,4) = [ 0.000 ; 0.500 ; 0.500 ] 
    atomchg(5) = atomselect(1); atompos(1:3,5) = [ 0.500 ; 0.000 ; 0.000 ] 
    atomchg(6) = atomselect(1); atompos(1:3,6) = [ 0.500 ; 0.000 ; 0.500 ] 
    atomchg(7) = atomselect(1); atompos(1:3,7) = [ 0.500 ; 0.500 ; 0.000 ] 
    atomchg(8) = atomselect(1); atompos(1:3,8) = [ 0.500 ; 0.500 ; 0.500 ]     
end

if isequal(inStructure,'Ah_t')
    cryst = 'cub'
    Nnu = 10
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.250 ; 0.250 ; 0.250 ] 
    atomchg(3) = atomselect(1); atompos(1:3,3) = [ 0.250 ; 0.250 ; 0.750 ] 
    atomchg(4) = atomselect(1); atompos(1:3,4) = [ 0.250 ; 0.750 ; 0.250 ] 
    atomchg(5) = atomselect(1); atompos(1:3,5) = [ 0.250 ; 0.750 ; 0.750 ]  
    atomchg(6) = atomselect(1); atompos(1:3,6) = [ 0.750 ; 0.250 ; 0.250 ] 
    atomchg(7) = atomselect(1); atompos(1:3,7) = [ 0.750 ; 0.250 ; 0.750 ] 
    atomchg(8) = atomselect(1); atompos(1:3,8) = [ 0.750 ; 0.750 ; 0.250 ] 
    atomchg(9) = atomselect(1); atompos(1:3,9) = [ 0.750 ; 0.750 ; 0.750 ]   
    atomchg(10) = atomselect(1); atompos(1:3,10) = [ 0.500 ; 0.500 ; 0.500 ]     
end

if isequal(inStructure,'A_p')
    cryst = 'cub'
    Nnu = 8
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.000 ; 0.667 ; 0.000 ] 
    atomchg(3) = atomselect(1); atompos(1:3,3) = [ 0.500 ; 0.000 ; 0.000 ] 
    atomchg(4) = atomselect(1); atompos(1:3,4) = [ 0.500 ; 0.667 ; 0.000 ] 
    atomchg(5) = atomselect(1); atompos(1:3,5) = [ 0.250 ; 0.500 ; 0.000 ]  
    atomchg(6) = atomselect(1); atompos(1:3,6) = [ 0.250 ; 0.167 ; 0.000 ] 
    atomchg(7) = atomselect(1); atompos(1:3,7) = [ 0.750 ; 0.500 ; 0.000 ] 
    atomchg(8) = atomselect(1); atompos(1:3,8) = [ 0.750 ; 0.167 ; 0.000 ]     
end


if isequal(inStructure,'A_c')
    cryst = 'cub'
    Nnu = 4
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.500 ; 0.000 ; 0.000 ] 
    atomchg(3) = atomselect(1); atompos(1:3,3) = [ 0.250 ; 0.250 ; 0.000 ] 
    atomchg(4) = atomselect(1); atompos(1:3,4) = [ 0.750 ; 0.250 ; 0.000 ]     
end


if isequal(inStructure,'A1')
    cryst = 'fcc'
    Nnu = 1; nunit = 1
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    filldensity = 74    
end

if isequal(inStructure,'Ah_Gas')
    cryst = 'cub'
    Nnu = 2
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.500 ; 0.500 ; 0.500 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.510 ; 0.510 ; 0.510 ]     
end

if isequal(inStructure,'A1_sc')
    cryst = 'cub'
    Nnu = 4
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.000 ; 0.500 ; 0.500 ] 
    atomchg(3) = atomselect(1); atompos(1:3,3) = [ 0.500 ; 0.000 ; 0.500 ] 
    atomchg(4) = atomselect(1); atompos(1:3,4) = [ 0.500 ; 0.500 ; 0.000 ] 
end

if isequal(inStructure,'A2')
    cryst = 'bcc'
    Nnu = 1
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    filldensity = 68     
end

if isequal(inStructure,'A3')
    cryst = 'hex'
    cbya = sqrt(8/3)
    Nnu = 2
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ]  
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.667 ; 0.333 ; 0.500 ]  
    filldensity = 74    
end

if isequal(inStructure,'A4')
    cryst = 'fcc'
    Nnu = 2
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.250 ; 0.250 ; 0.250 ]     
    filldensity = 34    
end


if isequal(inStructure,'B1')
    cryst = 'fcc'
    Nnu = 2
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(2); atompos(1:3,2) = [ 0.500 ; 0.500 ; 0.500 ]     
    filldensity = 52    
end

if isequal(inStructure,'B2')
    cryst = 'cub'
    Nnu = 2
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(2); atompos(1:3,2) = [ 0.500 ; 0.500 ; 0.500 ] 
    filldensity = 68     
end

if isequal(inStructure,'B3')
    cryst = 'fcc'
    Nnu = 2
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(2); atompos(1:3,2) = [ 0.250 ; 0.250 ; 0.250 ] 
    filldensity = 34      
end

if isequal(inStructure,'B4')
    cryst = 'hex'
    cbya = sqrt(8/3)
    Nnu = 4
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.667 ; 0.333 ; 0.500 ]  
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.333 ; 0.667 ; 0.000 ] 
    atomchg(3) = atomselect(2); atompos(1:3,3) = [ 0.667 ; 0.333 ; 0.875 ]  
    atomchg(4) = atomselect(2); atompos(1:3,4) = [ 0.333 ; 0.667 ; 0.375 ]  
    filldensity = 34     
end


if isequal(inStructure,'B_Z')
    cryst = 'fcc'
    Nnu = 4
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.500 ; 0.500 ; 0.500 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.750 ; 0.750 ; 0.750 ]     
    atomchg(3) = atomselect(2); atompos(1:3,3) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(4) = atomselect(2); atompos(1:3,4) = [ 0.250 ; 0.250 ; 0.250 ]     
end

if isequal(inStructure,'C1')
    cryst = 'fcc'
    Nnu = 3
    atomchg(1)  = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2)  = atomselect(2); atompos(1:3,2) = [ 0.250 ; 0.250 ; 0.250 ]     
    atomchg(3)  = atomselect(2); atompos(1:3,3) = [ 0.750 ; 0.750 ; 0.750 ]      
end

if isequal(inStructure,'C2')
    cryst = 'fcc'
    Nnu = 2
    atomchg(1)  = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2)  = atomselect(2); atompos(1:3,2) = [ 0.380 ; 0.380 ; 0.380 ]            
end

if isequal(inStructure,'C3')
    cryst = 'cub'
    Nnu = 6
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.250 ; 0.250 ; 0.250 ]
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.250 ; 0.750 ; 0.750 ]
    atomchg(3) = atomselect(1); atompos(1:3,3) = [ 0.750 ; 0.250 ; 0.750 ]
    atomchg(4) = atomselect(1); atompos(1:3,4) = [ 0.750 ; 0.750 ; 0.250 ] 
    atomchg(5) = atomselect(2); atompos(1:3,5) = [ 0.000 ; 0.000 ; 0.000 ]  
    atomchg(6) = atomselect(2); atompos(1:3,6) = [ 0.500 ; 0.500 ; 0.500 ]  
    filldensity = 26     
end

if isequal(inStructure,'C9')
    cryst = 'fcc'
    Nnu = 6
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ]  
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.250 ; 0.250 ; 0.250 ] 
    atomchg(3) = atomselect(2); atompos(1:3,3) = [ 0.125 ; 0.125 ; 0.125 ]
    atomchg(4) = atomselect(2); atompos(1:3,4) = [ 0.375 ; 0.375 ; 0.125 ]
    atomchg(5) = atomselect(2); atompos(1:3,5) = [ 0.375 ; 0.125 ; 0.375 ]
    atomchg(6) = atomselect(2); atompos(1:3,6) = [ 0.125 ; 0.375 ; 0.375 ]   
end

if isequal(inStructure,'C15')
    cryst = 'fcc'
    Nnu = 6
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ]  
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.250 ; 0.250 ; 0.250 ]                     
    atomchg(3) = atomselect(2); atompos(1:3,3) = [ 0.625 ; 0.625 ; 0.625 ]
    atomchg(4) = atomselect(2); atompos(1:3,4) = [ 0.625 ; 0.875 ; 0.875 ]
    atomchg(5) = atomselect(2); atompos(1:3,5) = [ 0.875 ; 0.625 ; 0.875 ]
    atomchg(6) = atomselect(2); atompos(1:3,6) = [ 0.875 ; 0.875 ; 0.625 ]  
    filldensity = 56      
end

if isequal(inStructure,'E2_1')
    cryst = 'cub'
    Nnu = 5
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ]  
    atomchg(2) = atomselect(2); atompos(1:3,2) = [ 0.500 ; 0.500 ; 0.500 ]  
    atomchg(3) = atomselect(3); atompos(1:3,3) = [ 0.000 ; 0.500 ; 0.500 ] 
    atomchg(4) = atomselect(3); atompos(1:3,4) = [ 0.500 ; 0.000 ; 0.500 ] 
    atomchg(5) = atomselect(3); atompos(1:3,5) = [ 0.500 ; 0.500 ; 0.000 ]
end

if isequal(inStructure,'E_C')
    cryst = 'cub'
    Nnu = 6
    atomchg(1) = atomselect(3); atompos(1:3,1) = [ 0.250 ; 0.250 ; 0.250 ]
    atomchg(2) = atomselect(3); atompos(1:3,2) = [ 0.250 ; 0.750 ; 0.750 ]
    atomchg(3) = atomselect(3); atompos(1:3,3) = [ 0.750 ; 0.250 ; 0.750 ]
    atomchg(4) = atomselect(3); atompos(1:3,4) = [ 0.750 ; 0.750 ; 0.250 ] 
    atomchg(5) = atomselect(1); atompos(1:3,5) = [ 0.000 ; 0.000 ; 0.000 ]  
    atomchg(6) = atomselect(2); atompos(1:3,6) = [ 0.500 ; 0.500 ; 0.500 ]   
end


if isequal(inStructure,'E2_X')
    cryst = 'cub'; Nnu  =  12;
    atomchg(1)  = atomselect(1); atompos(1:3,1)  = [ 0.000 ; 0.000 ; 0.000 ]
    atomchg(2)  = atomselect(1); atompos(1:3,2)  = [ 0.500 ; 0.500 ; 0.500 ]
    atomchg(3)  = atomselect(2); atompos(1:3,3)  = [ 0.500 ; 0.000 ; 0.250 ]
    atomchg(4)  = atomselect(2); atompos(1:3,4)  = [ 0.000 ; 0.500 ; 0.750 ]
    atomchg(5)  = atomselect(3); atompos(1:3,5)  = [ 0.250 ; 0.250 ; 0.125 ]
    atomchg(6)  = atomselect(3); atompos(1:3,6)  = [ 0.750 ; 0.750 ; 0.125 ] 
    atomchg(7)  = atomselect(3); atompos(1:3,7)  = [ 0.250 ; 0.750 ; 0.875 ]  
    atomchg(8)  = atomselect(3); atompos(1:3,8)  = [ 0.750 ; 0.250 ; 0.875 ] 
    atomchg(9)  = atomselect(3); atompos(1:3,9)  = [ 0.750 ; 0.750 ; 0.625 ]
    atomchg(10) = atomselect(3); atompos(1:3,10) = [ 0.250 ; 0.250 ; 0.625 ] 
    atomchg(11) = atomselect(3); atompos(1:3,11) = [ 0.750 ; 0.250 ; 0.375 ]  
    atomchg(12) = atomselect(3); atompos(1:3,12) = [ 0.250 ; 0.750 ; 0.375 ]  
end

if isequal(inStructure,'L1_2')
    cryst = 'cub'
    Nnu = 4; 
    atomchg(1) = atomselect(2); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ]   
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.000 ; 0.500 ; 0.500 ] 
    atomchg(3) = atomselect(1); atompos(1:3,3) = [ 0.500 ; 0.000 ; 0.500 ] 
    atomchg(4) = atomselect(1); atompos(1:3,4) = [ 0.500 ; 0.500 ; 0.000 ]
end


if isequal(inStructure,'Q_B3')
    cryst = 'cub'
    Nnu = 8;
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ] 
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.000 ; 0.500 ; 0.500 ] 
    atomchg(3) = atomselect(2); atompos(1:3,3) = [ 0.500 ; 0.000 ; 0.500 ] 
    atomchg(4) = atomselect(2); atompos(1:3,4) = [ 0.500 ; 0.500 ; 0.000 ] 
    atomchg(5) = atomselect(3); atompos(1:3,5) = [ 0.250 ; 0.250 ; 0.250 ]
    atomchg(6) = atomselect(3); atompos(1:3,6) = [ 0.250 ; 0.750 ; 0.750 ]
    atomchg(7) = atomselect(4); atompos(1:3,7) = [ 0.750 ; 0.250 ; 0.750 ]
    atomchg(8) = atomselect(4); atompos(1:3,8) = [ 0.750 ; 0.750 ; 0.250 ]
    filldensity = 34
end

if isequal(inStructure,'T_C15')
    cryst = 'fcc'
    Nnu = 6;
    atomchg(1) = atomselect(1); atompos(1:3,1) = [ 0.000 ; 0.000 ; 0.000 ]  
    atomchg(2) = atomselect(1); atompos(1:3,2) = [ 0.250 ; 0.250 ; 0.250 ]                     
    atomchg(3) = atomselect(2); atompos(1:3,3) = [ 0.625 ; 0.625 ; 0.625 ]
    atomchg(4) = atomselect(3); atompos(1:3,4) = [ 0.625 ; 0.875 ; 0.875 ]
    atomchg(5) = atomselect(2); atompos(1:3,5) = [ 0.875 ; 0.625 ; 0.875 ]
    atomchg(6) = atomselect(3); atompos(1:3,6) = [ 0.875 ; 0.875 ; 0.625 ]  
end

quiet silent


quiet normal
disp '*** Calculate the atomic energies and charges' 
quiet silent

%for (i=1:87)
%    espdofz(i,1:8) = [-0.410;  -0.122; -0.0702; ...]
%end

% atomic energies
    espdofz( 1,1:8) = [-0.233;  0.000;  0.000;  1; 0; 0; 031; 2.0] %120
    espdofz( 2,1:8) = [-0.570;  0.000;  0.000;  2; 0; 0; 050; 2.0]     
    espdofz( 3,1:8) = [-0.106;  0.000;  0.000;  1; 0; 0; 152; 1.7]    
    espdofz( 4,1:8) = [-0.206;  0.000;  0.000;  2; 0; 0; 111; 1.6] 
    espdofz( 5,1:8) = [-0.344; -0.137;  0.000;  2; 1; 0; 084; 1.6]    
    espdofz( 6,1:8) = [-0.500; -0.199;  0.000;  2; 2; 0; 076; 1.7]     
    espdofz( 7,1:8) = [-0.676; -0.266;  0.000;  2; 3; 0; 071; 2.0] 
    espdofz( 8,1:8) = [-0.871; -0.338;  0.000;  2; 4; 0; 066; 2.0]      
    espdofz( 9,1:8) = [-1.087; -0.416;  0.000;  2; 5; 0; 057; 1.5]     
    espdofz(10,1:8) = [-1.323; -0.498;  0.000;  2; 6; 0; 050; 2.0]   %058
    espdofz(11,1:8) = [-0.103;  0.000;  0.000;  1; 0; 0; 186; 1.9]   %166 %186
    espdofz(12,1:8) = [-0.175;  0.000;  0.000;  2; 0; 0; 160; 2.0]     
    espdofz(13,1:8) = [-0.287; -0.103;  0.000;  2; 1; 0; 143; 1.85] %1.85
    espdofz(14,1:8) = [-0.398; -0.153;  0.000;  2; 2; 0; 111; 1.9] 
    espdofz(15,1:8) = [-0.512; -0.206;  0.000;  2; 3; 0; 107; 2.0]      
    espdofz(16,1:8) = [-0.631; -0.261;  0.000;  2; 4; 0; 105; 2.0] 
    espdofz(17,1:8) = [-0.754; -0.320;  0.000;  2; 5; 0; 075; 1.4] %102   %075  
    espdofz(18,1:8) = [-0.883; -0.382;  0.000;  2; 6; 0; 106; 2.0]  
    espdofz(19,1:8) = [-0.089;  0.000;  0.000;  1; 0; 0; 227; 2.1]      
    espdofz(20,1:8) = [-0.141;  0.000;  0.000;  2; 0; 0; 197; 2.1]  
    espdofz(21,1:8) = [-0.156;  0.000; -0.131;  2; 0; 1; 161; 2.3]   
    espdofz(22,1:8) = [-0.167;  0.000; -0.170;  2; 0; 2; 145; 2.2]     
    espdofz(23,1:8) = [-0.176;  0.000; -0.204;  2; 0; 3; 131; 2.2] 
    espdofz(24,1:8) = [-0.150;  0.000; -0.118;  2; 0; 4; 124; 2.3]  
    espdofz(25,1:8) = [-0.191;  0.000; -0.267;  2; 0; 5; 137; 2.2]      
    espdofz(26,1:8) = [-0.198;  0.000; -0.296;  2; 0; 6; 124; 2.1]  
    espdofz(27,1:8) = [-0.204;  0.000; -0.322;  2; 0; 7; 125; 2.0] 
    espdofz(28,1:8) = [-0.210;  0.000; -0.349;  2; 0; 8; 125; 2.0]     
    espdofz(29,1:8) = [-0.172;  0.000; -0.202;  1; 0;10; 128; 2.0]   
    espdofz(30,1:8) = [-0.222;  0.000; -0.399;  2; 0;10; 133; 1.8]     
    espdofz(31,1:8) = [-0.328; -0.102;  0.000;  2; 1; 0; 122; 2.0] 
    espdofz(32,1:8) = [-0.427; -0.150;  0.000;  2; 2; 0; 122; 2.1]  %120  2.1
    espdofz(33,1:8) = [-0.523; -0.197;  0.000;  2; 3; 0; 119; 2.0]   
    espdofz(34,1:8) = [-0.621; -0.246;  0.000;  2; 4; 0; 120; 2.0]      
    espdofz(35,1:8) = [-0.720; -0.295;  0.000;  2; 5; 0; 120; 2.0] 
    espdofz(37,1:8) = [-0.085;  0.000;  0.000;  1; 0; 0; 248; 3.3]     
    espdofz(38,1:8) = [-0.132;  0.000;  0.000;  2; 0; 0; 215; 2.6] 
    espdofz(39,1:8) = [-0.151;  0.000; -0.109;  2; 0; 1; 178; 2.3] 
    espdofz(40,1:8) = [-0.162;  0.000; -0.151;  2; 0; 2; 159; 2.3]     
    espdofz(41,1:8) = [-0.144;  0.000; -0.125;  2; 0; 3; 143; 2.7]
    espdofz(42,1:8) = [-0.148;  0.000; -0.153;  2; 0; 4; 136; 2.6]    
    espdofz(47,1:8) = [-0.157;  0.000; -0.299;  1; 0;10; 144; 2.9]  %2.8 %144
    espdofz(49,1:8) = [-0.290; -0.102;  0.000;  2; 1; 0; 163; 3.0]  % 3 geschätzt       
    espdofz(50,1:8) = [-0.369; -0.144;  0.000;  2; 2; 0; 151; 3.3]     
    espdofz(51,1:8) = [-0.446; -0.186;  0.000;  2; 3; 0; 139; 2.0]    
    espdofz(52,1:8) = [-0.521; -0.227;  0.000;  2; 4; 0; 138; 2.0]      
    espdofz(53,1:8) = [-0.596; -0.268;  0.000;  2; 5; 0; 139; 2.0]  
    espdofz(55,1:8) = [-0.008;  0.000;  0.000;  1; 0; 0; 244; 2.0]   
    espdofz(56,1:8) = [-0.119;  0.000;  0.000;  2; 0; 0; 217; 3.4]     
    espdofz(57,1:8) = [-0.132;  0.000; -0.141;  2; 0; 1; 187; 2.0]      
    espdofz(58,1:8) = [-0.134;  0.000; -0.140;  2; 0; 1; 182; 2.0]  % Ce Start 4f
    espdofz(59,1:8) = [-0.124;  0.000; -0.001;  2; 0; 0; 182; 2.2] 
    espdofz(60,1:8) = [-0.126;  0.000; -0.001;  2; 0; 0; 181; 2.2] 
    espdofz(61,1:8) = [-0.127;  0.000; -0.001;  2; 0; 0; 181; 2.3] 
    espdofz(62,1:8) = [-0.128;  0.000; -0.001;  2; 0; 0; 180; 2.3] 
    espdofz(63,1:8) = [-0.129;  0.000; -0.001;  2; 0; 0; 193; 2.0]     
    espdofz(64,1:8) = [-0.143;  0.000; -0.127;  2; 0; 1; 179; 2.3]   
    espdofz(65,1:8) = [-0.132;  0.000; -0.001;  2; 0; 0; 177; 2.3] 
    espdofz(66,1:8) = [-0.133;  0.000; -0.001;  2; 0; 0; 175; 2.3] 
    espdofz(67,1:8) = [-0.134;  0.000; -0.001;  2; 0; 0; 174; 2.3] 
    espdofz(68,1:8) = [-0.135;  0.000; -0.001;  2; 0; 0; 173; 2.3] 
    espdofz(69,1:8) = [-0.136;  0.000; -0.001;  2; 0; 0; 172; 2.3] 
    espdofz(70,1:8) = [-0.137;  0.000; -0.001;  2; 0; 0; 194; 2.7] 
    espdofz(71,1:8) = [-0.155;  0.000; -0.103;  2; 0; 1; 172; 2.3]   % Lu Ende 4f
    espdofz(73,1:8) = [-0.175;  0.000; -0.182;  2; 0; 3; 143; 2.0]  
    espdofz(74,1:8) = [-0.181;  0.000; -0.220;  2; 0; 4; 137; 3.1]   %162
    espdofz(79,1:8) = [-0.162;  0.000; -0.305;  1; 0;10; 144; 2.5]      
    espdofz(80,1:8) = [-0.205;  0.000; -0.453;  2; 0;10; 150; 2.0]    
    espdofz(81,1:8) = [-0.285; -0.101;  0.000;  2; 1; 0; 170; 2.0]    
    espdofz(82,1:8) = [-0.357; -0.142;  0.000;  2; 2; 0; 175; 2.0]  
    espdofz(84,1:8) = [-0.494; -0.218;  0.000;  2; 4; 0; 168; 2.3]       
    espdofz(87,1:8) = [ 0.000;  0.000;  0.000;  0; 0; 0; 100; 2.0]  
    espdofz(88,1:8) = [ 0.000;  0.000;  0.000; -1; 0; 0; 100; 2.0] 
    espdofz(89,1:8) = [ 0.000;  0.000;  0.000;  1; 0; 0; 100; 2.0]     
    eatom = 0
    nexp = 0 
    lmax(1:Nnu) = 2
    for (iat=1:Nnu)
         e_at(iat,1) = espdofz(atomchg(iat),1)
         e_at(iat,2) = espdofz(atomchg(iat),2)
         e_at(iat,3) = espdofz(atomchg(iat),3)
         nval(iat)   = sum(espdofz(atomchg(iat),4:6))
         qvalnl(iat,1) = espdofz(atomchg(iat),4)
         qvalnl(iat,2) = espdofz(atomchg(iat),5)
         qvalnl(iat,3) = espdofz(atomchg(iat),6)    
         ratom(iat) = espdofz(atomchg(iat),7)
         eatom = eatom + sum(e_at(iat,1:3).*espdofz(atomchg(iat),4:6))      
         if (abs(e_at(iat,3)) > 0)
             lmax(iat) = 3
         end             
         nexp = nexp + espdofz(atomchg(iat),8)^2
    end 
    nexp = -sqrt(nexp/Nnu)
%    nexp = -Nnu/sum(1./espdofz(atomchg(:),8))
%    nexp = -mean(espdofz(atomchg(:),7).^3.*espdofz(atomchg(:),8)/mean(espdofz(atomchg(:),7).^3))
    nvalence = sum(nval)

    
qs=sum(qvalnl(:,1))
qp=sum(qvalnl(:,2))
qd=sum(qvalnl(:,3))
%nest = (2*(qs+qp) + 3.5*(qs+qp)*qd + 5*qd)/(qs+qp+(qs+qp)*qd + qd)
nest = -(2*(qs+qp) + 0.01*3.5*(qs+qp)*qd)/(qs+qp+0.01*(qs+qp)*qd)

quiet normal
disp '*** Building the real and reciprocal lattice and atomic positions' 
quiet silent
    

kboltz = 8.617e-5   % Boltzmannkonstante in eV/K
rgas = 8.314 % universelle gaskonstante in J/(mol*K)
ehart = 0.03675 % 1/Hartree-Energie
anull = 0.529 % Bohr-Radius
gpa = 160.2  % Umrechnung ev/Angström^3 auf GPa
kgprol = 1.6606
mpros = 1000
evnm = 1240
evinkjmol = 96.485
blin = 135 % emp. Faktor für Lindemannsche Schmelzpktformel

GPoint = [0.0, 0.0, 0.0]  % G
if (isequal(cryst,'fcc'))
   t1 = [0.0, 0.5, 0.5]
   t2 = [0.5, 0.0, 0.5]
   t3 = [0.5, 0.5, 0.0]
   APoint = [0.50, 0.50, 0.50]  % L
   BPoint = [0.00, 1.00, 0.00]  % X
   CPoint = [0.50, 0.50, 0.00]  % K 
   b1titel = ['Bandindex L - G - X - K - Projected DOS']
end
if (isequal(cryst,'cub'))
   t1 = [1.0, 0.0, 0.0]
   t2 = [0.0, 1.0, 0.0]
   t3 = [0.0, 0.0, 1.0]
   APoint = [0.5, 0.5, 0.0]  % M
   BPoint = [0.5, 0.0, 0.0]  % X
   CPoint = [0.5, 0.5, 0.5]  % R
   b1titel = ['Bandindex M - G - X - R - Projected DOS']
end
if (isequal(cryst,'bcc'))
   t1 = [-0.5,  0.5,  0.5]
   t2 = [ 0.5, -0.5,  0.5]
   t3 = [ 0.5,  0.5, -0.5]
   APoint = [0.25, 0.25, 0.25] % P
   BPoint = [0.00, 0.50, 0.00] % H
   CPoint = [0.25, 0.25, 0.00] % N   
   b1titel = ['Bandindex P - G - H - N - Projected DOS']
end
if (isequal(cryst,'hex'))
   t1 = [ 0.5, -sqrt(3)/2, 0.0]
   t2 = [ 0.5,  sqrt(3)/2, 0.0]
   t3 = [ 0.0,  0.0, cbya]
   APoint = [0.50, 0.50, 0.00] % K
   BPoint = [0.00, 0.00, 0.50] % A  
   CPoint = [0.50, 0.00, 0.50] % L    
   b1titel = ['Bandindex K - G - A - L - Projected DOS']
end
if (isequal(cryst,'ort'))
   t1 = [1.0, 0.0, 0.0]
   t2 = [0.0, bbya, 0.0]
   t3 = [0.0, 0.0, cbya]
   APoint = [0.5, 0.5, 0.0]
   BPoint = [0.5, 0.5, 0.5]   
   CPoint = [0.0, 0.5, 0.0]
   b1titel = ['Bandindex S - G - R - Y - Projected DOS']
end
if (isequal(cryst,'tet'))
   t1 = [1.0, 0.0, 0.0]
   t2 = [0.0, 1.0, 0.0]
   t3 = [0.0, 0.0, cbya]
   APoint = [0.5, 0.5, 0.0]
   BPoint = [0.5, 0.5, 0.5]   
   CPoint = [0.0, 0.5, 0.0]
   b1titel = ['Bandindex M - G - A - X - Projected DOS']
end

kpoint = [APoint;GPoint;BPoint;CPoint]
nkband = 4
%
% determine the real spacelattice vectors
T = [t1;t2;t3]; volume = det(T)
t23 = [t2(2)*t3(3)-t2(3)*t3(2),t2(3)*t3(1)-t2(1)*t3(3),t2(1)*t3(2)-t2(2)*t3(1)]
t31 = [t3(2)*t1(3)-t3(3)*t1(2),t3(3)*t1(1)-t3(1)*t1(3),t3(1)*t1(2)-t3(2)*t1(1)]
t12 = [t1(2)*t2(3)-t1(3)*t2(2),t1(3)*t2(1)-t1(1)*t2(3),t1(1)*t2(2)-t1(2)*t2(1)]
% Determine the set of reciprocal space and k-Vectors 
b1 = t23./volume; b2 = t31./volume ; b3 = t12./volume
B = [b1;b2;b3]
for (i=1:nkband)
    kpoint(i,:) = [sum(kpoint(i,1).*b1),sum(kpoint(i,2).*b2),sum(kpoint(i,3).*b3)]
end
%
% build a super cell and determine the vectors between the neighbored cells
ncut = 1; iq = 0
for (h=-ncut:ncut); for (k=-ncut:ncut); for (l=-ncut:ncut)
    iq = iq + 1
    vecpos(1:3,iq) = sum([h.*t1;k.*t2;l.*t3])        
end; ;end; end; nq = iq 
%
% determine the real atomic positions within the unit cell and the super cell
xatomcor(1:3,1:Nnu,1:nq) = 0
for (iat=1:Nnu) 
    atomcor(:,iat) = sum([atompos(1,iat).*t1;atompos(2,iat).*t2;atompos(3,iat).*t3]) 
 if (isequal(cryst,'fcc'))   
    atomcor = atompos
 end
end 
for (iat=1:Nnu)
   for (iq = 1:nq)
      xatomcor(:,iat,iq) = vecpos(:,iq) + atomcor(:,iat) 
   end
end
%
% Determine the distances between the atoms
dmin(1:Nnu,1:Nnu) = 10000; dist(1:Nnu,1:Nnu,1:nq) = 0;
dx(1:Nnu,1:Nnu,1:nq) = 0; dy(1:Nnu,1:Nnu,1:nq) = 0; dz(1:Nnu,1:Nnu,1:nq) = 0
for (iat = 1:Nnu); for (jat = 1:Nnu); for (iq=1:nq)
   dx(iat,jat,iq) = atomcor(1:1,iat) - xatomcor(1:1,jat,iq) 
   dy(iat,jat,iq) = atomcor(2:2,iat) - xatomcor(2:2,jat,iq) 
   dz(iat,jat,iq) = atomcor(3:3,iat) - xatomcor(3:3,jat,iq) 
   dist(iat,jat,iq) = sqrt(dx(iat,jat,iq)^2 + dy(iat,jat,iq)^2 + dz(iat,jat,iq)^2)
   if (dist(iat,jat,iq) > 0) 
      if (dmin(iat,jat) > dist(iat,jat,iq))
          dmin(iat,jat) = dist(iat,jat,iq)
      end        
   end 
end; end; end 
% determine the nearest neighbor
sr(1:Nnu,1:Nnu,1:nq) = 0; near(1:Nnu) = 0
for (iat = 1:Nnu);for (jat = 1:Nnu); for (iq=1:nq)
    if (dist(iat,jat,iq) > 0)        
      sr(iat,jat,iq) = (min(dmin(iat,:))/dist(iat,jat,iq))^12  
    end      
    near(iat) = near(iat) + sr(iat,jat,iq)
end; end; end

%
% Muffin-Tin Näherung
if (isequal(alat,0))
    for (iat=1:Nnu)
        for (jat=1:Nnu)
             aest(iat,jat)=(ratom(iat)+ratom(jat))/dmin(iat,jat)
        end
    end
    alat = max(max(aest))/100
end               
alat = alat/anull
%
% tightbinding parameters

quiet normal
disp '*** Calculate the interatomic matrix elements ' 
quiet silent

pdr(1:3,1:3) = 0; Vsig(1:3,1:3) = 0; Vpi(1:3,1:3)= 0; Ylm(1:3,1:5) = 0; 
Elmg00(1:Nnu,1:9,1:Nnu,1:9,1:nkband) = 0; psk(1:9,1:9) = 0


pdr(1:2,1:2) = 2; pdr(3,1:2) = 3.5; pdr(1:2,3) = 3.5; pdr(3,3) = 5    
nss = -1.32; nsp = 1.42; nxx =  2.22; nxy = -0.63 
% theoretical values
% nss = -9*pi^2/64; nsp = 9*pi^2/32*sqrt(1-16/(3*pi^2)); nxx=21*pi^2/64; nxy = -3*pi^2/32
ndds = -35.6; nddp = 19.2; nsd = -4.68; npds = -4.37; npdp = 2.02

Vsig(1,1) =  nss/alat^pdr(1,1)
Vsig(1,2) =  nsp/alat^pdr(1,2)
Vsig(2,1) = -nsp/alat^pdr(2,1)
Vsig(2,2) =  nxx/alat^pdr(2,2)
Vsig(1,3) =  nsd/alat^pdr(1,3)
Vsig(3,1) = -nsd/alat^pdr(3,1)
Vsig(3,2) =  npds/alat^pdr(3,2)
Vsig(2,3) = -npds/alat^pdr(2,3)
Vsig(3,3) =  ndds/alat^pdr(3,3)
Vpi0(2,2) =  nxy/alat^pdr(2,2)
Vpi0(2,3) =  npdp/alat^pdr(2,3)
Vpi0(3,2) =  npdp/alat^pdr(3,2)  
Vpi0(3,3) =  nddp/alat^pdr(3,3)

for (ik=1:nkband)
    kv = [kpoint(ik,1),kpoint(ik,2),kpoint(ik,3)]  
 %       kv = [kpoint(1,1),kpoint(1,2),kpoint(1,3)]  
    kxp = 2*pi*kv(1); kyp = 2*pi*kv(2); kzp = 2*pi*kv(3)  
    for (iat = 1:Nnu)
        for (jat = 1:Nnu)
            for (iq=1:nq)   
               if (sr(iat,jat,iq) > 0.01)            
                    if (dist(iat,jat,iq) > 0) 
                        g00r = cos((dx(iat,jat,iq)*kxp+dy(iat,jat,iq)*kyp+dz(iat,jat,iq)*kzp))                    
                        g00i = j*sin((dx(iat,jat,iq)*kxp+dy(iat,jat,iq)*kyp+dz(iat,jat,iq)*kzp))*0 
                        g00 = g00r + g00i
                        dl = dx(iat,jat,iq)/dist(iat,jat,iq)
                        dm = dy(iat,jat,iq)/dist(iat,jat,iq)
                        dn = dz(iat,jat,iq)/dist(iat,jat,iq) 
                        Ylm(1,1) = 1                        % s
                        Ylm(2,1) = dl                       % px
                        Ylm(2,2) = dm                       % py
                        Ylm(2,3) = dn                       % pz
                        Ylm(3,1) = sqrt(3)/2*(dl*dl-dm*dm)  % dx2-y2
                        Ylm(3,2) = dn*dn-(dm*dm+dl*dl)/2    % d3z2-r2
                        Ylm(3,3) = sqrt(3)*dl*dm            % dxy
                        Ylm(3,4) = sqrt(3)*dl*dn            % dxz
                        Ylm(3,5) = sqrt(3)*dm*dn            % dzy      
         if (1<0)
                        psk(2,2) = 1-dl*dl                            % px - px
                        psk(2,3) = -dl*dm                             % px - py
                        psk(2,4) = -dl*dn                             % px - pz
                        psk(3,3) = 1-dm*dm                            % py - py
                        psk(3,4) = -dm*dn                             % py - pz                        
                        psk(4,4) = 1-dn*dn                            % pz - pz
                        psk(2,5) = dl*(1-dl*dl+dm*dm)                 % px - dx2-y2
                        psk(2,6) = -sqrt(3)*dl*dn*dn                  % px - d3z2-r2                     
                        psk(2,7) = dm*(1-2*dl*dl)                     % px - dxy
                        psk(2,8) = dn*(1-2*dl*dl)                     % px - dxz
                        psk(2,9) = -2*dl*dm*dn                        % px - dzy
                        psk(3,5) = dm*(1-dm*dm+dl*dl)                 % py - dx2-y2
                        psk(3,6) = -sqrt(3)*dm*dn*dn                  % py - d3z2-r2                     
                        psk(3,7) = dl*(1-2*dm*dm)                     % py - dxy
                        psk(3,8) = -2*dl*dm*dn                        % py - dxz
                        psk(3,9) = dn*(1-2*dm*dm)                     % py - dzy 
                        psk(4,5) = -dn*(dl*dl-dm*dm)                  % pz - dx2-y2
                        psk(4,6) = sqrt(3)*dn*(dl*dl+dm*dm)           % pz - d3z2-r2                     
                        psk(4,7) = -2*dl*dm*dn                        % pz - dxy
                        psk(4,8) = dl*(1-2*dn*dn)                     % pz - dxz
                        psk(4,9) = dm*(1-2*dn*dn)                     % pz - dzy    
                        psk(5,5) = dl*dl+dm*dm-(dl*dl-dm*dm)^2        % dx2-y2 - dx2-y2
                        psk(5,6) = dn*dn*(dm*dm-dl*dl)^2              % dx2-y2 - d3z2-r2                     
                        psk(5,7) = 2*dl*dm*(dm*dm-dl*dl)              % dx2-y2 - dxy
                        psk(5,8) = dn*dl*(1-2*(dl*dl-dm*dm))          % dx2-y2 - dxz
                        psk(5,9) = dn*dm*(1+2*(dl*dl-dm*dm))          % dx2-y2 - dzy    
                        psk(6,6) = 3*dn*dn*(dm*dm+dl*dl)              % d3z2-r2 - d3z2-r2                     
                        psk(6,7) = 2*dl*dm*dn*dn                      % d3z2-r2 - dxy
                        psk(6,8) = dl*dn*(dm*dm+dl*dl-dn*dn)          % d3z2-r2 - dxz
                        psk(6,9) = dm*dn*(dm*dm+dl*dl-dn*dn)          % d3z2-r2 - dzy    
                        psk(7,7) = dl*dl+dm*dm-4*dl*dl*dm*dm          % dxy - dxy
                        psk(7,8) = dm*dn*(1-4*dl*dl)                  % dxy - dxz
                        psk(7,9) = dl*dn*(1-4*dm*dm)                  % dxy - dzy
                        psk(8,8) = dl*dl+dn*dn-4*dl*dl*dn*dn          % dxz - dxz
                        psk(8,9) = dl*dm*(1-4*dn*dn)                  % dxz - dzy
                        psk(9,9) = dm*dm+dn*dn-4*dm*dm*dn*dn          % dzy - dzy                        
                        for (is=1:9)
                            for (it=is+1:9)
                                 psk(it,is) = psk(is,it)
                            end
                        end
end                      
                        ilv = 0                      
                        for (il=1:lmax(iat))
                            ilm = 2*il - 1
                            for (im=1:ilm)   
                                ilv = ilv + 1
                                jlv = 0                
                                for (jl=1:lmax(jat))
                                    jlm = 2*jl - 1
                                    for (jm=1:jlm)
                                        jlv = jlv + 1
                                        pv = isequal(il,jl)*isequal(im,jm)
                                        dr = 1/dist(iat,jat,iq)^pdr(il,jl)
                                        
Esig = dr*sr(iat,jat,iq)*Ylm(il,im)*Ylm(jl,jm)*Vsig(il,jl)*g00*(4/near(iat))^0.5
Epi0 = dr*sr(iat,jat,iq)*(pv - Ylm(il,im)*Ylm(jl,jm))*Vpi0(il,jl)*g00 
%if (ik == 1)
%    if (il == 2)
%        if (im == 3)
%            if (jl == 2)
%                if (jm == 1)
%                    quiet normal
%                    test = [Epi0, dr,sr(iat,jat,iq),pv,Ylm(il,im),Ylm(jl,jm),Vpi0(il,jl),g00]
%                    quiet silent
%                end
%            end
%        end
%    end
%end 
        
                                        
                                        Elmg00(iat,ilv,jat,jlv,ik) = Elmg00(iat,ilv,jat,jlv,ik) + dr*sr(iat,jat,iq)*Ylm(il,im)*Ylm(jl,jm)*Vsig(il,jl)*g00*(4/near(iat))^0.5
%                                        Elmg00(iat,ilv,jat,jlv,ik) = Elmg00(iat,ilv,jat,jlv,ik) + dr*sr(iat,jat,iq)*psk(ilv,jlv)*Vpi0(il,jl)*g00 %*(4/near(iat))^0.5                                      
                                        Elmg00(iat,ilv,jat,jlv,ik) = Elmg00(iat,ilv,jat,jlv,ik) + dr*sr(iat,jat,iq)*(pv - Ylm(il,im)*Ylm(jl,jm))*Vpi0(il,jl)*g00 %*(4/near(iat))^0.5                                     
                                    end
                                end
                            end
                        end                      
                    end    
                end
            end                 
       end
    end
end

quiet normal
disp '*** Calculate the band structure eigenvalues for each k ' 
quiet silent
    
ndim = sum(lmax(:).*lmax(:))    
igo = 0
bandek(1:ndim,1:nkband) = 0
ekksum(1:nkband*ndim) = 0
for (ik=1:nkband)     
    H = zeros(ndim)  
    ihmat = 0
    for (iat = 1:Nnu)
        ilv = 0
        for (il=1:lmax(iat))
            ilm = il+il-1
            for (im = 1:ilm)
                ilv = ilv + 1
                ihmat = ihmat + 1
                jhmat = 0              
                for (jat = 1:Nnu)
                  jlv = 0
                  for (jl = 1:lmax(jat))
                    jlm  = jl + jl - 1
                    for (jm = 1:jlm)
                       jlv = jlv + 1
                       jhmat = jhmat + 1

       H(ihmat,jhmat) = e_at(iat,il)*isequal(iat,jat)*isequal(il,jl)*isequal(im,jm) + Elmg00(iat,ilv,jat,jlv,ik)

                end
              end
            end
        end
      end
    end     
    
    eev = sort(real(eig(H)/ehart)) 
    for (g=1:ndim)
        igo = igo + 1
        bandek(g,ik) = eev(g)  
        ekksum(igo)  = eev(g)        
    end
    if (isequal(ik,2))  
        nfermi=floor(nvalence/2)
        if (isequal(nfermi,0))
            nfermi = 1
        end 
        if (isequal(nfermi,ndim))
            nfermi = ndim - 1
        end
        gapdirect = eev(nfermi+1) - eev(nfermi) 
    end
end


nkbpan = round((nkbfit-1)/3)
nkbfit = 3*nkbpan + 1
kp1= linspace(0,1,nkbpan+1)
bndfit(1:ndim,1:nkbfit) = 0
for (g=1:ndim)
    for (ipan=1:3)
        db = bandek(g,ipan+1) - bandek(g,ipan)
        bndfit(g,(ipan-1)*nkbpan+1:ipan*nkbpan+1) = bandek(g,ipan)+db*(0.5-cos(pi*kp1)/2)
    end
end


quiet normal
disp '*** Calculate the density of states' 
quiet silent


% DOS linear interpolation between Symmetry-Points
jgo=0
for (g=1:ndim)
    jgo = jgo + 1
    ekkfit(jgo) = bandek(g,2)
end
nmix = 10000/ndim/7
for (ig=1:nmix)
    kpo = (ig/nmix)^0.5
    for (g=1:ndim)
        db1 = bandek(g,1) - bandek(g,2)
        db2 = bandek(g,3) - bandek(g,2)
        db3 = bandek(g,4) - bandek(g,2)
        ekf(1) = bandek(g,2)+db1*kpo %(0.5-cos(pi*kpo)/2)
        ekf(2) = bandek(g,2)+db2*kpo %(0.5-cos(pi*kpo)/2)
        ekf(3) = bandek(g,2)+db3*kpo %(0.5-cos(pi*kpo)/2)
        ekf(4) = (ekf(1) + ekf(2) + ekf(3))/3
        ekf(5) = (ekf(2) + ekf(3))/2
        ekf(6) = (ekf(1) + ekf(3))/2
        ekf(7) = (ekf(1) + ekf(2))/2        
        for (ia=1:7)
            jgo = jgo + 1
            ekkfit(jgo) = ekf(ia)
        end
    end
end


ngo = ndim*nkband
mgo = jgo
efsort = sort(ekkfit)
nfermi = floor(mgo*nvalence/ndim/2)
if (isequal(nfermi,0))
   nfermi = 1
end 
if (isequal(nfermi,mgo))
   nfermi = mgo - 1
end 
efermi = (efsort(nfermi) + efsort(nfermi+1))/2
bandgap = efsort(nfermi+1)-efsort(nfermi)
bandwidth = efsort(1) - efermi + bandgap/2
%
% Set the Fermi-Energy to zero
bandek = bandek - efermi
ekksum = ekksum - efermi
efsort = efsort - efermi
bndfit = bndfit - efermi
emaxx = max(efsort)
emins = min(efsort)
eminy = floor(emins/5)*5
emaxy = ceil(emaxx/5)*5

nmeshd = 200 
defbox = (emaxy-eminy)/(nmeshd-1)  
efspace=linspace(eminy,emaxy,nmeshd) 
dosfit = hist(efsort,efspace) 
normdos = sum(dosfit)
dosfit = dosfit./normdos*ndim*2/defbox
eweight=1./(exp((efspace)/(kboltz*1))+1)


for (imesh=1:nmeshd-1)
    if (efspace(imesh) < 0)
        if (efspace(imesh+1) > 0)
           dosatef = (dosfit(imesh)+dosfit(imesh+1))/2
           nfermi = imesh
           dosmov = dosfit(imesh+1)*dosfit(imesh)
        end
    end 
    if (efspace(imesh) == 0)       
        dosatef = dosfit(imesh)
        nfermi = imesh   
        dosmov = dosfit(imesh+1)*dosfit(imesh)     
    end    
end

gapratio = bandgap/gapdirect
substtype = ' metal'
gaptype = 'no'
if (isequal(dosatef,0))
    substtype = 'n insulator'
    if (bandgap < 2.5)
       substtype = ' semiconductor'
    end 
    gaptype = 'a direct'
    if (gapratio < 1)
        gaptype = 'an indirect'
    end
end

quiet normal
disp '*** Adding f-electrons to the density-of-states if necessary' 
quiet silent
fourf = [1,3,4,5,6,7,7,9,10,11,12,13,14,14]
for iat=1:Nnu
    if (atomchg(iat) > 57)
        if (atomchg(iat) < 72)
            eoff = -0.3  % f-energy level
            fsharp = 10  % sharpness of the f-band
            efocc = eoff/ehart - efermi
            focc = fourf(atomchg(iat) - 57)
            fuoc = 14 - focc
            dosplf = focc*sqrt(fsharp*focc/pi)*exp(-fsharp*focc*(efspace-efocc).^2) + fuoc*sqrt(fsharp*fuoc/pi)*exp(-fsharp*fuoc*(efspace+1/2*efocc).^2)
            dosfit = dosfit + dosplf
            eatom = eatom + eoff*focc
        end 
    end
end


quiet normal
disp '*** Calculate the bonding properties and physical properties' 
quiet silent

%quiet normal
alat = alat*anull
molmasse = sum(0.0063.*atomchg.*atomchg + 2.0886.*atomchg -  1.3453) % sum(2.006*atomchg)
volnull = volume*alat^3
dichte = molmasse/volnull*kgprol
eatom = eatom/ehart
chargeint = sum(qvalnl(:,1)  + qvalnl(:,2).*0.33) %+ qvalnl(:,3).*(10-qvalnl(:,3))/12)
rhoint = chargeint/volnull
evalence = sum((efspace+efermi).*dosfit.*eweight)*defbox

Ebond = (evalence-eatom) % + 11*volnull*(chargeint/volnull)^(5/3)
rwigner = (3*volnull/4/pi/Nnu)^(1/3)
bulkmodul = -gpa*Ebond/volnull*(nexp+1)^2
betab = bulkmodul/gpa
shearmodul = 0.4*bulkmodul
c11 = bulkmodul+4/3*shearmodul
youngmodul = 9*bulkmodul*shearmodul/(3*bulkmodul+shearmodul)
debye = 177.5*(Nnu)^(1/3)*(volnull)^(1/6)*sqrt(bulkmodul/molmasse)
vl=sqrt((3*bulkmodul+4*shearmodul)/(3*dichte))*mpros
vt=sqrt(shearmodul/dichte)*mpros
schall = (1/3*(2/vt^3+1/vl^3))^(-1/3)
Enull = -betab*volnull/(nexp+1)^2/Nnu
ediss = -Enull 
rhoint=(bulkmodul/1900)^(3/5)
grueneisen= -1/6-nexp 
mohs = 10*tanh(2*bulkmodul/500)+0.2
vdiss =exp(-1/nexp+log(volnull))
tschmelz = 3e-5*debye^2*(volnull/Nnu)^(2/3)*molmasse/Nnu

%quiet silent

tspace = linspace(10,600,100)
meanmol = molmasse/Nnu
meanvol = volnull/Nnu

for it=(1:100)   
    xdt = linspace(debye/tspace(it)/1000,debye/tspace(it),100); dtx=xdt(2)-xdt(1) 
    cvt(it) = 9*rgas*(tspace(it)/debye)^3*sum(xdt.^4.*exp(xdt)./(exp(xdt)-1).^2)*dtx
    if (tspace(it) < 293)
        irt = it
    end
    ffermi=1./(exp((efspace)/(kboltz*tspace(it)))+1)
    chgmov(it) =  sum(dosfit.*abs(ffermi-eweight))*defbox + 1e-30 %+ dosatef*defbox % %rhoint*sum(dosfit.*(1-ffermi).*ffermi)*defbox/nvalence 
    tauex = 5
    phono(it) = ((tspace(it)/debye)^tauex*sum(xdt.^tauex./((exp(xdt)-1).*(1-exp(-xdt))))*dtx)
end

cvt = cvt/meanmol*1000
cvrt = cvt(irt)
alphal = grueneisen/3.*cvt(irt)/bulkmodul*dichte
% rho = (1/(1e-30 + chgmov) + 20*phono + 2*sum(qvalnl(:,3).*(10-qvalnl(:,3))))*2.0e-9
%rho = 5e-8/dosatef + 5e-8*phono./rhoint % + 2*sum(qvalnl(:,3).*(10-qvalnl(:,3)))*6.0e-9

% rho = 1e-8 + 1e-8*phono./rhoint

v=linspace(volnull*0.9,volnull*1.1,50)
etotal = betab*volnull/(nexp+1).*(v/volnull).^(nexp+1).*(log(v./volnull)-1/(nexp+1))
e1tot =  betab*volnull/(nexp+1).*(v/volnull).^(nexp+2)
e2tot = -betab*volnull/(nexp+1).*(v/volnull).^(nexp+1)
e3tot = -betab*volnull/(nexp+1)^2.*(v/volnull).^(nexp+1)
ptotal = -betab*(v/volnull).^(nexp).*log(v./volnull)*gpa
ptotal1 = (betab*(v/volnull).^(nexp))*gpa
ptotal2 = (-betab*(v/volnull).^(nexp+1))*gpa

quiet normal
disp '*** Calculate the optical properties' 
quiet silent
% base of the optical properties ist the (Pseudo-)Joint-DOS 
pjdos(1:nmeshd) = 0
for (jmesh = 1:nfermi)
    for (gmesh = nfermi+1:nmeshd)
        pjdos(gmesh-jmesh) = pjdos(gmesh-jmesh) + dosfit(jmesh)*dosfit(gmesh)*defbox*defbox
    end
end   
% using a constant oscillator strength we calculated the transition probability eps2
% and via Kronig-Kramer relationship the dielectric/suceptibility functions
fosz=22000/((ndim - nvalence/2)*(nvalence/2))
hnuf = linspace(defbox,emaxy-eminy+defbox,nmeshd)
omega = hnuf
hmin = round(0.8/defbox)
hmax = round(12.4/defbox)
if (hmax > nmeshd)
    hmax = nmeshd
end
eps2 = pi/volnull*fosz*pjdos./omega.^2
dshift = defbox/2
for (iw = 1:nmeshd)
    omegas = omega(iw)+dshift
    chi(iw) = 2/pi*sum(omega.*eps2./(omega.^2-omegas^2))
end
eps1 = 1 + chi
% which give us the frquency dependend optical indices, reflection and transmission
nindex = sqrt(0.5*(sqrt(eps1(1)^2+eps2(1)^2)+eps1(1)))
kindex = sqrt(0.5*(sqrt(eps1(1)^2+eps2(1)^2)-eps1(1)))
reflection = ((nindex-1)/(nindex+1))^2
transmission = 1 - reflection
normeps2 = eps2(hmin:hmax)/max(eps2(hmin:hmax))
normeps1 = eps1(hmin:hmax)/max(abs(eps1(hmin:hmax)))
lambda = evnm/hnuf(hmin:hmax)   
dosfit(1) = 0; dosfit(nmeshd) = 0 
%
% Calculate the optical gap using the Joint-DOS - could be very different in case of f-electrons
optigap = defbox
for (iw=2:nmeshd-1)
    if (pjdos(nmeshd-iw)>0.0001)
        optigap = hnuf(nmeshd-iw)
    end
end
absorpedge = evnm/optigap
if (isequal(optigap,defbox))
    optigap = 0
end

% elect. conductivity resp resistance
% Elektron-Phonon interaction increases with T and increases R
rhoep = 10e-8*phono
% Electron-electron: Large portion of movable electrons decreases R 
rhoee = 1.e-9*Nnu/(dosatef*defbox+chgmov) 
% s-d electron scattering: Increases R with filled-to-unfilled d-states and collison with movable electrons
rhosd = 4e-8*sum(qvalnl(:,3).*(10-qvalnl(:,3))).*chgmov/Nnu
rho = (rhoep + rhoee + rhosd)
alrhot = (rho(irt+1)-rho(irt))/(rho(irt)*(tspace(irt+1)-tspace(irt)))
qrese = 2.0e-8*1/rho.*tspace 
%wfunc = tanh(0.54*tspace/debye)
%whigh = 1./tspace
%wlow  =  exp(debye./(2*tspace)
qresp = 6e-6*schall*dichte*cvt.*(1*0+(debye/tspace))
kappa = (qrese + qresp)

quiet normal
disp '*** Plot the results' 
quiet silent


subplot(2,5,[1,2,6,7])    
knull = linspace(nkbpan+1,nkbpan+1,2)
kzwo = linspace(2*nkbpan+1,2*nkbpan+1,2)
kend = linspace(3*nkbpan+1,3*nkbpan+1,2)
enull = linspace(eminy,emaxy,2)
elevl = linspace(0,0,2)
klevl = linspace(1,nkbfit,2)
mxds = max(dosfit)
dosplt = dosfit.*nkbpan/mxds
plot(1:nkbfit,bndfit(:,1:nkbfit),'-b',knull,enull,'k',kzwo,enull,'k',klevl,elevl,'k',kend,enull,'k',nkbfit+dosplt,efspace,'b')
xlim([1,4*nkbpan+1])
xlabel(b1titel)
ylabel('Energy in eV')
ylim([eminy,emaxy])
drawnow
    
subplot(2,5,[3:4])
semilogx (lambda,normeps2,lambda,normeps1)
ylabel(' eps1, eps2')
xlabel('Photon wavelength in nm')
xlim([100,1400])
drawnow

subplot(2,5,5)
if (optigap > 0)
    plot (tspace,1/rho*1e-6)
    ylabel('sigma(T) in MS/m')
end
if (isequal(optigap,0))
    plot (tspace,rho*1e8)
    ylabel('rho(T) in *10^-8 Ohm/m')
end
xlabel('T in K')
drawnow

subplot(2,5,8)
plot (v/Nnu,etotal*evinkjmol/Nnu)
ylabel('E in kJ/mol per atom')
xlabel('Volume in A^3 per atom')
drawnow

subplot(2,5,9)
plot (tspace,cvt)
ylabel('cv(T) in J/(kg*K)')
xlabel('T in K')
drawnow

subplot(2,5,10)
symb = ['ovxxodooovxxxdooovxx**********xdooovxx**********xdooovxx************************************']
for (iat=1:Nnu)
    for (iq=1:nq)
    isincell = sqrt((xatomcor(1,iat,iq)-0.5)^2+(xatomcor(2,iat,iq)-0.5)^2+(xatomcor(3,iat,iq)-0.5)^2)  
    if (isincell < sqrt(1))
       plot3(xatomcor(1,iat,iq),xatomcor(2,iat,iq),xatomcor(3,iat,iq),symb(atomchg(iat)),'linewidth',1)
       hold on 
%       drawnow    
    end
    end   
end
for (iat=1:Nnu)   
    for (jat=1:Nnu)
        for (iq=1:nq)
            if (sr(iat,jat,iq) == 1)
                bondx = linspace(atomcor(1,iat),xatomcor(1,jat,iq),2)
                bondy = linspace(atomcor(2,iat),xatomcor(2,jat,iq),2)  
                bondz = linspace(atomcor(3,iat),xatomcor(3,jat,iq),2)                   
                plot3(bondx,bondy,bondz,'linewidth',1)  
                hold on
                drawnow
            end    
        end    
    end   
end
xlim([0,1])
ylim([0,1])
zlim([0,1])
hold off
view(15,20)
hold off
ylabel('a')
xlabel('b')
zlabel('c')
drawnow


print('Plot_result.png')

quiet normal

printf ('Results of the calculation:\n')
printf ('============================\n')
printf ('The structuretype is %s with a %s crystal structure and a lattice constant of %.3f A.\n', inStructure, cryst, alat)
printf ('The equilibrium volume is %.2f A^3 with %i atoms per unit cell. \n', volnull, Nnu)
printf ('The molare mass is %.2f g/mol with a density of %.3f kg/l. \n', molmasse, dichte)
printf ('The bulk modulus is %.1f GPa, the shear modulus is %.1f GPa and the Young modulus %.1f GPa \n', bulkmodul, shearmodul, youngmodul)
printf ('The linear thermal expansion coefficient is %.2f µm/(m*K) and the Grüneisen factor is %.2f.\n', alphal, grueneisen)
printf ('The Mohs hardness is %.0f and the Debye temperature is %.1f K. \n', mohs, debye)
printf ('The average speed of sound is %.1f m/s (longitudinal velocity: %.1f m/s, transversal velocity: %.1f m/s). \n', schall, vl, vt)
printf ('The specific heat capacity at room temperature is %.1f J/(kg*K). \n', cvrt)
printf ('The melting temperature is about %.1f K. \n', tschmelz)
printf ('The vaporization enthalpie is %.1f kJ/mol. \n', ediss*evinkjmol )
printf ('The specific thermal conductivity at room temerature is %.1f W/mK. \n', kappa(irt))
printf ('The optical absorption edge lies at %.1f nm. The refractive index is %.2f.  \n' , absorpedge, nindex)
printf ('The degree of reflexion is %.1f Prozent. The degree of transmission is %.1f Prozent.\n', reflection*100, transmission*100 )
printf ('The Fermi enegry is %.2f eV. The valence band has a broadening of %.2f eV. \n', efermi, bandwidth)
printf ('The minimum bandgap of the host is %.2f eV. The optical band gap ist %.2f. \n',  bandgap, optigap)
printf ('The charge density at the Fermi level is %.2f e/eV. \n',  dosatef)
printf ('The specific electrical resistance is %.2e Ohm/m resp. of %.2e S/m at room temperature \n',  rho(irt), 1/rho(irt))
printf ('with a temperature coefficient of %.2e per K at room temperature \n',  alrhot)
printf ('The substance is a%s with %s bandgap \n', substtype, gaptype)
if (dosatef > 0)
printf ('and a interstitial electron charge density of %.1f*10-2 e/A^3. \n', rhoint*100)  
end


