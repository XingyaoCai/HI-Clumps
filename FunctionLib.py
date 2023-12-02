import illustris_python as il
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
basePath ='./sims.TNG/TNG50-1/output'

def calculate_R90(Rmin, Rmax, Masses,Coordinates,CenterOfMass):
    precision = 1
    max_iter = 1000


    left, right = Rmin, Rmax

    for i in range(max_iter):

        mid = (left + right) / 2


        ratio =  mass_ratio_within_radius(mid, Masses, Coordinates, CenterOfMass)


        if ratio < 0.9:
            left = mid

        else:
            right = mid


        if abs(right - left) < precision or i == max_iter - 1:
            R90 = (left + right) / 2
            break

    return R90

def mass_ratio_within_radius(Radius, Masses, Coordinates, CenterOfMass):

    dist_to_com = np.sqrt(np.sum((Coordinates - CenterOfMass)**2, axis=1))


    within_radius_idx = np.where(dist_to_com <= Radius)[0]


    total_mass_within_radius = np.sum(Masses[within_radius_idx])

    total_mass = np.sum(Masses)

    mass_ratio = total_mass_within_radius / total_mass

    return mass_ratio

def Binary_DataIter(Data_Set1,Data_Set2):
    for i in range(0,min(len(Data_Set1),len(Data_Set2))):
        yield Data_Set1[i],Data_Set2[i]


def drawsubhalo(i,sp,savepath,basePath,String):
    gas_data=il.snapshot.loadSubhalo(basePath, sp, i, 'gas', fields=['Coordinates','NeutralHydrogenAbundance','Masses'])
    if len(gas_data)==1:return -1

    GasCoordinates = gas_data['Coordinates']
    GasAbundance = gas_data['NeutralHydrogenAbundance']
    Masses = gas_data['Masses']
    NeutralHydrogenAbundance=Masses*GasAbundance



    GasCoorMin = GasCoordinates.min(axis=0)
    GasCoorMax = GasCoordinates.max(axis=0)
    CircleCenter = (GasCoorMin + GasCoorMax) / 2
    plt.figure(figsize=(12, 12), dpi=200)

    plt.style.use("dark_background")
    bars=(GasCoorMax-GasCoorMin)
    h, _, _, image = plt.hist2d(GasCoordinates[:, 1],
                            GasCoordinates[:, 2],
                            weights=NeutralHydrogenAbundance,

                            norm=mpl.colors.LogNorm(),
                            bins=(bars[1],bars[2]))
    plt.xlim([GasCoorMin[1],GasCoorMax[1]])
    plt.ylim([GasCoorMin[2],GasCoorMax[2]])
    plt.xlabel('y [ckpc/h]')
    plt.ylabel('z [ckpc/h]')
    plt.text(0,0,String)
    #draw_circle = plt.Circle(( CircleCenter[1],  CircleCenter[2]), 150, fill=False)
    #plt.gcf().gca().add_artist(draw_circle)
    plt.colorbar(image)
    plt.gca().set_aspect(1)
    plt.savefig(savepath+'/{}.png'.format(i))


def subhaloID2subfindID (subhaloID,sublink):
    subhaloID_array=sublink['SubhaloID']
    subfindIndex=np.where(subhaloID_array==subhaloID)
    subfindID=sublink['SubfindID'][subfindIndex]

    return subfindID

def DescendantSubfindID(subfindID,snapshot_num,sublink):
    tree_fields=['SubfindID','SnapNum','DescendantID']
    subfindid_snapnum_tree=il.sublink.loadTree(basePath, snapshot_num,subfindID, fields=tree_fields, onlyMPB=True)
    if subfindid_snapnum_tree['DescendantID'][0]==-1:return -1
    descendant_subfind_id=subhaloID2subfindID(subfindid_snapnum_tree['DescendantID'][0],sublink)[0]
    return descendant_subfind_id



def Total_Spin_Cal(Spin_Array):
    total_spin=0
    for spin in Spin_Array:
        total_spin+=(spin)**2
    total_spin=total_spin**(1/2)
    return total_spin


def mass_within_radius(Radius, Masses, Coordinates, CenterOfMass):
    # 计算每个恒星到质心的距离
    dist_to_com = np.sqrt(np.sum((Coordinates - CenterOfMass)**2, axis=1))

    # 找到距离质心小于R_200的所有恒星的下标
    within_radius_idx = np.where(dist_to_com <= Radius)[0]

    # 计算在R_200范围内的恒星总质量
    total_mass_within_radius = np.sum(Masses[within_radius_idx])

    return total_mass_within_radius


def Cold_Gas_Mass(snap_num,subhalo_id):
    gas_fields=['ElectronAbundance','InternalEnergy','Masses']
    subhalo_data=il.snapshot.loadSubhalo(basePath,snap_num,subhalo_id,0,fields=gas_fields)
    if len(subhalo_data)==1:
        return 0

    x_e=subhalo_data['ElectronAbundance'].astype(np.float64)
    internal_energy=subhalo_data['InternalEnergy'].astype(np.float64)
    gas_cell_masses=subhalo_data['Masses'].astype(np.float64)

    m_p=1.673E-24
    X_H=0.76
    unit_switching=1E10
    mean_molecular_weight=4*m_p/(1+3*X_H+4*X_H*x_e)
    k_B=1.38E-16
    gas_cell_temperature_in_Kelvin=2/3*internal_energy/k_B*unit_switching*mean_molecular_weight
    cold_gas_mass=0
    gas_mass=gas_cell_masses.sum()

    for temperature,gas_cell_mass in Binary_DataIter(gas_cell_temperature_in_Kelvin,gas_cell_masses):
        if temperature<10000:cold_gas_mass+=gas_cell_mass

    subhalo_mass=subhalos_data[subhalo_id]
    return cold_gas_mass, subhalo_mass

def drawhalo(i,sp,savepath,basePath):
    gas_data=il.snapshot.loadHalo(basePath, sp, i, 'gas', fields=['Coordinates','NeutralHydrogenAbundance','Masses'])

    GasCoordinates = gas_data['Coordinates']
    GasAbundance = gas_data['NeutralHydrogenAbundance']
    Masses = gas_data['Masses']
    NeutralHydrogenAbundance=Masses*GasAbundance

    GasCoorMin = GasCoordinates.min(axis=0)
    GasCoorMax = GasCoordinates.max(axis=0)
    CircleCenter = (GasCoorMin + GasCoorMax) / 2
    plt.figure(figsize=(12, 12), dpi=200)

    plt.style.use("dark_background")
    bars=(GasCoorMax-GasCoorMin)
    h, _, _, image = plt.hist2d(GasCoordinates[:, 1],
                            GasCoordinates[:, 2],
                            weights=NeutralHydrogenAbundance,

                            norm=mpl.colors.LogNorm(),
                            bins=(bars[1],bars[2]))
    plt.xlim([GasCoorMin[1],GasCoorMax[1]])
    plt.ylim([GasCoorMin[2],GasCoorMax[2]])
    plt.xlabel('y [ckpc/h]')
    plt.ylabel('z [ckpc/h]')
    #draw_circle = plt.Circle(( CircleCenter[1],  CircleCenter[2]), 150, fill=False)
    #plt.gcf().gca().add_artist(draw_circle)
    plt.colorbar(image)
    plt.gca().set_aspect(1)
    plt.savefig(savepath+'/{}.png'.format(i))




