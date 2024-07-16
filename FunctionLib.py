import illustris_python as il
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


class Cell:
    def __init__(self):
        self.Coordinate = np.array([0.0, 0.0, 0.0])
        self.Mass = 0.0


class Cloud:
    def __init__(self, x_range, y_range):
        x_min, x_max = x_range
        y_min, y_max = y_range
        x_centers = np.arange(x_min + 0.0625, x_max, 0.125)
        y_centers = np.arange(y_min + 0.0625, y_max, 0.125)
        xx, yy = np.meshgrid(x_centers, y_centers)
        self.cells = [Cell() for _ in range(xx.size)]
        self.x_range = x_range
        self.y_range = y_range
        self.assign_coordinates(xx.flatten(), yy.flatten())

    def assign_coordinates(self, x_coords, y_coords):
        for i, cell in enumerate(self.cells):
            cell.Coordinate = np.array([0.0, y_coords[i], x_coords[i]])

    def assign_mass(self, coordinates: np.ndarray, masses: np.ndarray):
        x_min, x_max = self.x_range
        y_min, y_max = self.y_range
        x_bars = int((x_max - x_min) / 0.125)
        for i in range(coordinates.shape[0]):
            coord = coordinates[i]
            mass = masses[i]
            x_index = int((coord[0] - x_min) / 0.125)
            y_index = int((coord[1] - y_min) / 0.125)

            index = x_index * x_bars + y_index
            if (index >= len(self.cells)):
                continue
            if (index < 0):
                continue
            self.cells[index].Mass += mass

    def visualize(self, savePath):
        plt.figure(figsize=(20, 20), dpi=300)
        plt.style.use("dark_background")
        x_range, y_range = self.x_range, self.y_range
        plt.hist2d([cell.Coordinate[1] for cell in self.cells], [cell.Coordinate[2] for cell in self.cells], norm=mpl.colors.LogNorm(
        ), cmin=0.00005, bins=[int((y_range[1]-y_range[0])*2), int((x_range[1]-x_range[0])*2)], weights=[cell.Mass for cell in self.cells])
        plt.xlabel('y [ckpc/h]')
        plt.ylabel('z [ckpc/h]')
        plt.ylim(self.x_range[0], self.x_range[1])
        plt.xlim(self.y_range[0], self.y_range[1])
        plt.axes().get_xaxis().set_visible(False)
        plt.axes().get_yaxis().set_visible(False)
        plt.axis('equal')
        plt.savefig(savePath)
        plt.close('all')

    def reset(self):
        for cell in self.cells:
            cell.Coordinate = np.array([0.0, 0.0, 0.0])
            cell.Mass = 0.0


def Binary_DataIter(Data_Set1, Data_Set2):
    for i in range(0, min(len(Data_Set1), len(Data_Set2))):
        yield Data_Set1[i], Data_Set2[i]


def draw(i, sp):
    GasCoordinates = il.snapshot.loadSubhalo(
        basePath, sp, i, 'gas', fields=['Coordinates'])
    # GasAbundance = il.snapshot.loadSubhalo(basePath, sp,i, 'gas', fields=['NeutralHydrogenAbundance'])
    Masses = il.snapshot.loadSubhalo(basePath, sp, i, 'gas', fields=['Masses'])
    # NeutralHydrogenAbundance=Masses*GasAbundance

    GasCoorMin = GasCoordinates.min(axis=0)
    GasCoorMax = GasCoordinates.max(axis=0)
    CircleCenter = (GasCoorMin + GasCoorMax) / 2
    plt.figure(figsize=(12, 12), dpi=200)

    plt.style.use("dark_background")
    bars = (GasCoorMax-GasCoorMin)
    h, _, _, image = plt.hist2d(GasCoordinates[:, 1],
                                GasCoordinates[:, 2],
                                weights=Masses,  # NeutralHydrogenAbundance,

                                norm=mpl.colors.LogNorm(),
                                bins=(bars[1], bars[2]))
    plt.xlim([GasCoorMin[1], GasCoorMax[1]])
    plt.ylim([GasCoorMin[2], GasCoorMax[2]])
    plt.xlabel('y [ckpc/h]')
    plt.ylabel('z [ckpc/h]')
    # draw_circle = plt.Circle(( CircleCenter[1],  CircleCenter[2]), 150, fill=False)
    # plt.gcf().gca().add_artist(draw_circle)
    plt.colorbar(image)
    plt.gca().set_aspect(1)
    plt.savefig('./img/pic-{}.png'.format(sp))
    plt.show()


def Draw_Circle(image_Path, Saved_Path, radius):
    image = Image.open(image_Path)
    draw = ImageDraw.Draw(image)
    width, height = image.size
    center_x, center_y = width // 2, height // 2
    pixel_radius = radius * convert_parameter
    circle_color = (255, 255, 255)
    border_width = 2  # 你可以根据需要调整边框宽度
    draw.ellipse(
        [
            (center_x - pixel_radius, center_y - pixel_radius),
            (center_x + pixel_radius, center_y + pixel_radius),
        ],
        outline=circle_color,
        width=border_width,
    )
    image.save(Saved_Path)
    image.close()


def baryon_mass_in_radius(Center, Radius, Coordinates, Masses):
    dist_to_com = np.sqrt(np.sum((Coordinates - Center)**2, axis=1))
    within_radius_idx = np.where(dist_to_com <= Radius)[0]
    total_mass_within_radius = np.sum(Masses[within_radius_idx])
    del dist_to_com, within_radius_idx
    return total_mass_within_radius


def dm_mass_in_radius(Center, Radius, Coordinates):
    dist_to_com = np.sqrt(np.sum((Coordinates - Center)**2, axis=1))
    within_radius_idx = np.where(dist_to_com <= Radius)[0]
    total_mass_within_radius = m_dm*within_radius_idx.shape[0]*h
    del dist_to_com, within_radius_idx
    return total_mass_within_radius


def total_mass_in_radius(gas_data_dict, dm_data_dict, stars_data_dict, bh_data_dict, Center, Radius):

    gas_mass_in = baryon_mass_in_radius(
        Center, Radius, gas_data_dict['Coordinates'], gas_data_dict['Masses'])
    dm_mass_in = dm_mass_in_radius(Center, Radius, dm_data_dict['Coordinates'])
    stars_mass_in = baryon_mass_in_radius(
        Center, Radius, stars_data_dict['Coordinates'], stars_data_dict['Masses'])
    bh_mass_in = baryon_mass_in_radius(
        Center, Radius, bh_data_dict['Coordinates'], bh_data_dict['Masses'])

    total_mass_in = gas_mass_in+dm_mass_in+stars_mass_in+bh_mass_in

    Return_Dict = dict()
    Return_Dict['Mass'] = total_mass_in
    Return_Dict['PartType0_Mass'] = gas_mass_in
    Return_Dict['PartType1_Mass'] = dm_mass_in
    Return_Dict['PartType4_Mass'] = stars_mass_in
    Return_Dict['PartType5_Mass'] = bh_mass_in
    return Return_Dict


def M_And_R_To_Vel(Mass, Radius):
    return np.sqrt(43007.1*Mass/Radius)


def V_Circ_At_Radius_From_Mass_Dict(Mass_Dict, Radius):
    Return_Dict = dict()
    Return_Dict['VCirc'] = M_And_R_To_Vel(Mass_Dict['Mass'], Radius)
    Return_Dict['PartType0_VCirc'] = M_And_R_To_Vel(
        Mass_Dict['PartType0_Mass'], Radius)
    Return_Dict['PartType1_VCirc'] = M_And_R_To_Vel(
        Mass_Dict['PartType1_Mass'], Radius)
    Return_Dict['PartType4_VCirc'] = M_And_R_To_Vel(
        Mass_Dict['PartType4_Mass'], Radius)
    Return_Dict['PartType5_VCirc'] = M_And_R_To_Vel(
        Mass_Dict['PartType5_Mass'], Radius)

    return Return_Dict


def SubhaloParticleRadius90(BasePath, Subhalo_Index, Snapshot_Number, BoxSize, PartType):

    particles_dict = il.snapshot.loadSubhalo(
        BasePath, Snapshot_Number, Subhalo_Index, PartType, fields=['Coordinates', 'Masses'])
    subhalo_center = il.groupcat.loadSingle(
        BasePath, Snapshot_Number, -1, Subhalo_Index)['SubhaloPos']
    Particles_Distance = batch_distance_calculation(
        particles_dict['Coordinates'], subhalo_center, BoxSize)
    particles_masses = particles_dict['Masses']
    Particles_Distance_Sorted = np.sort(Particles_Distance)
    Particles_Masses_Sorted = particles_masses[np.argsort(Particles_Distance)]
    Masses_Cumulative = np.cumsum(Particles_Masses_Sorted)
    Masses_Cumulative_Normalized = Masses_Cumulative/Masses_Cumulative[-1]
    Radius90 = Particles_Distance_Sorted[np.where(
        Masses_Cumulative_Normalized >= 0.9)[0][0]]
    return Radius90

def Binary_DataIter(Data_Set1,Data_Set2):
    for i in range(0,min(len(Data_Set1),len(Data_Set2))):
        yield Data_Set1[i],Data_Set2[i]


def calculate_coldgas_mass(snap_num, subhalo_id):
    gas_fields = ['ElectronAbundance', 'InternalEnergy', 'Masses']
    subhalo_data = il.snapshot.loadSubhalo(basePath, snap_num, subhalo_id, 0, fields=gas_fields)
    if len(subhalo_data) == 1:
        return 0

    x_e = subhalo_data['ElectronAbundance'].astype(np.float64)
    internal_energy = subhalo_data['InternalEnergy'].astype(np.float64)
    gas_cell_masses = subhalo_data['Masses'].astype(np.float64)

    m_p = 1.673E-24
    X_H = 0.76
    unit_switching = 1E10
    mean_molecular_weight = 4 * m_p / (1 + 3 * X_H + 4 * X_H * x_e)
    k_B = 1.38E-16
    gas_cell_temperature_in_Kelvin = 2 / 3 * internal_energy / k_B * unit_switching * mean_molecular_weight

    cold_gas_mask = np.where(gas_cell_temperature_in_Kelvin < 10000)
    cold_gas_mass = gas_cell_masses[cold_gas_mask].sum()
    gas_mass = gas_cell_masses.sum()

    return cold_gas_mass, gas_mass

def calculate_stars_mass(snap_num, subhalo_id):
    stars_fields = ['Masses']
    subhalo_data = il.snapshot.loadSubhalo(basePath, snap_num, subhalo_id, 4, fields=stars_fields)
    if len(subhalo_data) == 1:
        return 0



    stars_mass = subhalo_data.sum()

    return stars_mass

def Cold_Gas_Ratio(snap_num,subhalo_id):
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
    ratio=cold_gas_mass/gas_mass
    subhalo_mass=subhalos_data[subhalo_id]
    return ratio, subhalo_mass

def batch_distance_calculation(pos1, pos2, boxsize=None):
    # 计算每个点与目标点的距离
    dxs = np.abs(pos1[:, 0] - pos2[0])
    dys = np.abs(pos1[:, 1] - pos2[1])
    dzs = np.abs(pos1[:, 2] - pos2[2])

    # 考虑周期性边界条件
    if boxsize is not None:
        dxs2 = np.abs(boxsize - dxs)
        dys2 = np.abs(boxsize - dys)
        dzs2 = np.abs(boxsize - dzs)

        # 在周期性边界条件下，选择最短的距离
        dxs = np.minimum(dxs, dxs2)
        dys = np.minimum(dys, dys2)
        dzs = np.minimum(dzs, dzs2)

    # 计算三维空间点之间的欧氏距离
    distances = np.linalg.norm(np.vstack([dxs, dys, dzs]).T, axis=1)

    return distances


def compute_velocity_distance(PartType0_Dict, Subhalo_Dict, chunk_size=1000):
    N = PartType0_Dict['Coordinates'].shape[0]
    num_chunks = (N + chunk_size - 1) // chunk_size

    norm_velocity_list = []
    tang_velocity_list = []
    distance_list = []

    for i in tqdm(range(num_chunks)):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, N)

        chunk_coordinates = PartType0_Dict['Coordinates'][start:end]
        chunk_velocities = PartType0_Dict['Velocities'][start:end]

        relative_coordinates = chunk_coordinates - Subhalo_Dict['SubhaloPos']
        relative_distance = np.linalg.norm(relative_coordinates, axis=1)
        relative_velocity = chunk_velocities - Subhalo_Dict['SubhaloVel']

        norm_velocity = np.sum(
            relative_velocity * relative_coordinates, axis=1) / relative_distance
        tang_velocity = np.sqrt(
            np.sum(relative_velocity**2, axis=1) - norm_velocity**2)

        norm_velocity_list.append(norm_velocity)
        tang_velocity_list.append(tang_velocity)
        distance_list.append(relative_distance)

    norm_velocity_array = np.concatenate(norm_velocity_list)
    tang_velocity_array = np.concatenate(tang_velocity_list)
    distance_array = np.concatenate(distance_list)

    return norm_velocity_array, tang_velocity_array, distance_array


def plot_velocity_distance(distance, tang_velocity, num_bins=50, distance_range=(0, 700)):
    bin_edges = np.linspace(*distance_range, num=num_bins+1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    bin_indices = np.digitize(distance, bin_edges)

    mean_velocity = np.array(
        [tang_velocity[bin_indices == i].mean() for i in range(1, num_bins+1)])
    std_velocity = np.array([tang_velocity[bin_indices == i].std()
                            for i in range(1, num_bins+1)])

    plt.scatter(distance, tang_velocity, 0.5, 'grey', alpha=0.1)
    plt.errorbar(bin_centers, mean_velocity,
                 yerr=std_velocity, fmt='o-', color='red')

    plt.xlabel('Distance [kpc / h]')
    plt.ylabel('Tangential Velocity [km/s]')


def compute_velocity_distance_with_filter(PartType0_Dict, Subhalo_Dict, x_difference_threshold=50, chunk_size=1000):
    N = PartType0_Dict['Coordinates'].shape[0]
    num_chunks = (N + chunk_size - 1) // chunk_size

    norm_velocity_list = []
    tang_velocity_list = []
    distance_list = []

    for i in range(num_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, N)

        chunk_coordinates = PartType0_Dict['Coordinates'][start:end]
        chunk_velocities = PartType0_Dict['Velocities'][start:end]

        # 添加筛选条件：x 坐标与给定 Subhalo_Dict 中的 x 坐标相差在阈值以内
        x_condition = np.abs(
            chunk_coordinates[:, 0] - Subhalo_Dict['SubhaloPos'][0]) <= x_difference_threshold
        chunk_coordinates = chunk_coordinates[x_condition]
        chunk_velocities = chunk_velocities[x_condition]

        relative_coordinates = chunk_coordinates - Subhalo_Dict['SubhaloPos']
        relative_distance = np.linalg.norm(relative_coordinates, axis=1)
        relative_velocity = chunk_velocities - Subhalo_Dict['SubhaloVel']

        norm_velocity = np.sum(
            relative_velocity * relative_coordinates, axis=1) / relative_distance
        tang_velocity = np.sqrt(
            np.sum(relative_velocity**2, axis=1) - norm_velocity**2)

        norm_velocity_list.append(norm_velocity)
        tang_velocity_list.append(tang_velocity)
        distance_list.append(relative_distance)

    norm_velocity_array = np.concatenate(norm_velocity_list)
    tang_velocity_array = np.concatenate(tang_velocity_list)
    distance_array = np.concatenate(distance_list)

    return norm_velocity_array, tang_velocity_array, distance_array


def plot_velocity_distance(distance, tang_velocity, VCirc, VCirc_Gas, VCirc_Stars, VCirc_DM, num_bins=50, distance_range=(0, 700)):
    bin_edges = np.linspace(*distance_range, num=num_bins+1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    bin_indices = np.digitize(distance, bin_edges)

    mean_velocity = np.array(
        [tang_velocity[bin_indices == i].mean() for i in range(1, num_bins+1)])
    std_velocity = np.array([tang_velocity[bin_indices == i].std()
                            for i in range(1, num_bins+1)])

    plt.scatter(distance, tang_velocity, 0.5, 'grey', alpha=0.1)
    plt.errorbar(bin_centers, mean_velocity,
                 yerr=std_velocity, fmt='o-', color='red')

    r = np.arange(start=0, stop=500, step=10)
    plt.plot(r, VCirc, 'b', label='VCirc')
    plt.plot(r, VCirc_Gas, 'red', label='VCirc_Gas')
    plt.plot(r, VCirc_Stars, 'yellow', label='VCirc_Stars')
    plt.plot(r, VCirc_DM, 'green', label='VCirc_DM')

    plt.xlabel('r [kpc]')
    plt.ylabel('V_Circ [km / s]')
    plt.legend()
    plt.show()


def mass_ratio_within_radius(Radius, Masses, Coordinates, CenterOfMass):

    dist_to_com = np.sqrt(np.sum((Coordinates - CenterOfMass)**2, axis=1))

    within_radius_idx = np.where(dist_to_com <= Radius)[0]

    total_mass_within_radius = np.sum(Masses[within_radius_idx])

    total_mass = np.sum(Masses)

    mass_ratio = total_mass_within_radius / total_mass

    return mass_ratio


def Binary_DataIter(Data_Set1, Data_Set2):
    for i in range(0, min(len(Data_Set1), len(Data_Set2))):
        yield Data_Set1[i], Data_Set2[i]


def drawsubhalo(i, sp, savepath, basePath, String):
    gas_data = il.snapshot.loadSubhalo(basePath, sp, i, 'gas', fields=[
                                       'Coordinates', 'NeutralHydrogenAbundance', 'Masses'])
    if len(gas_data) == 1:
        return -1

    GasCoordinates = gas_data['Coordinates']
    GasAbundance = gas_data['NeutralHydrogenAbundance']
    Masses = gas_data['Masses']
    NeutralHydrogenAbundance = Masses*GasAbundance

    GasCoorMin = GasCoordinates.min(axis=0)
    GasCoorMax = GasCoordinates.max(axis=0)
    CircleCenter = (GasCoorMin + GasCoorMax) / 2
    plt.figure(figsize=(12, 12), dpi=200)

    plt.style.use("dark_background")
    bars = (GasCoorMax-GasCoorMin)
    h, _, _, image = plt.hist2d(GasCoordinates[:, 1],
                                GasCoordinates[:, 2],
                                weights=NeutralHydrogenAbundance,

                                norm=mpl.colors.LogNorm(),
                                bins=(bars[1], bars[2]))
    plt.xlim([GasCoorMin[1], GasCoorMax[1]])
    plt.ylim([GasCoorMin[2], GasCoorMax[2]])
    plt.xlabel('y [ckpc/h]')
    plt.ylabel('z [ckpc/h]')
    plt.text(0, 0, String)
    # draw_circle = plt.Circle(( CircleCenter[1],  CircleCenter[2]), 150, fill=False)
    # plt.gcf().gca().add_artist(draw_circle)
    plt.colorbar(image)
    plt.gca().set_aspect(1)
    plt.savefig(savepath+'/{}.png'.format(i))


def subhaloID2subfindID(subhaloID, sublink):
    subhaloID_array = sublink['SubhaloID']
    subfindIndex = np.where(subhaloID_array == subhaloID)
    subfindID = sublink['SubfindID'][subfindIndex]

    return subfindID


def DescendantSubfindID(subfindID, snapshot_num, sublink, basePath):
    tree_fields = ['SubfindID', 'SnapNum', 'DescendantID']
    subfindid_snapnum_tree = il.sublink.loadTree(
        basePath, snapshot_num, subfindID, fields=tree_fields, onlyMPB=True)
    if subfindid_snapnum_tree['DescendantID'][0] == -1:
        return -1
    descendant_subfind_id = subhaloID2subfindID(
        subfindid_snapnum_tree['DescendantID'][0], sublink)[0]
    return descendant_subfind_id


def Total_Spin_Cal(Spin_Array):
    total_spin = 0
    for spin in Spin_Array:
        total_spin += (spin)**2
    total_spin = total_spin**(1/2)
    return total_spin


def mass_within_radius(Radius, Masses, Coordinates, CenterOfMass):
    # 计算每个恒星到质心的距离
    dist_to_com = np.sqrt(np.sum((Coordinates - CenterOfMass)**2, axis=1))

    # 找到距离质心小于R_200的所有恒星的下标
    within_radius_idx = np.where(dist_to_com <= Radius)[0]

    # 计算在R_200范围内的恒星总质量
    total_mass_within_radius = np.sum(Masses[within_radius_idx])

    return total_mass_within_radius


def Cold_Gas_Mass(snap_num, subhalo_id, basePath, subhalos_data):
    gas_fields = ['ElectronAbundance', 'InternalEnergy', 'Masses']
    subhalo_data = il.snapshot.loadSubhalo(
        basePath, snap_num, subhalo_id, 0, fields=gas_fields)
    if len(subhalo_data) == 1:
        return 0

    x_e = subhalo_data['ElectronAbundance'].astype(np.float64)
    internal_energy = subhalo_data['InternalEnergy'].astype(np.float64)
    gas_cell_masses = subhalo_data['Masses'].astype(np.float64)

    m_p = 1.673E-24
    X_H = 0.76
    unit_switching = 1E10
    mean_molecular_weight = 4*m_p/(1+3*X_H+4*X_H*x_e)
    k_B = 1.38E-16
    gas_cell_temperature_in_Kelvin = 2/3*internal_energy / \
        k_B*unit_switching*mean_molecular_weight
    cold_gas_mass = 0
    gas_mass = gas_cell_masses.sum()

    for temperature, gas_cell_mass in Binary_DataIter(gas_cell_temperature_in_Kelvin, gas_cell_masses):
        if temperature < 10000:
            cold_gas_mass += gas_cell_mass

    subhalo_mass = subhalos_data[subhalo_id]
    return cold_gas_mass, subhalo_mass


def drawhalo(i, sp, savepath, basePath):
    gas_data = il.snapshot.loadHalo(basePath, sp, i, 'gas', fields=[
                                    'Coordinates', 'NeutralHydrogenAbundance', 'Masses'])

    GasCoordinates = gas_data['Coordinates']
    GasAbundance = gas_data['NeutralHydrogenAbundance']
    Masses = gas_data['Masses']
    NeutralHydrogenAbundance = Masses*GasAbundance

    GasCoorMin = GasCoordinates.min(axis=0)
    GasCoorMax = GasCoordinates.max(axis=0)
    CircleCenter = (GasCoorMin + GasCoorMax) / 2
    plt.figure(figsize=(12, 12), dpi=200)

    plt.style.use("dark_background")
    bars = (GasCoorMax-GasCoorMin)
    h, _, _, image = plt.hist2d(GasCoordinates[:, 1],
                                GasCoordinates[:, 2],
                                weights=NeutralHydrogenAbundance,

                                norm=mpl.colors.LogNorm(),
                                bins=(bars[1], bars[2]))
    plt.xlim([GasCoorMin[1], GasCoorMax[1]])
    plt.ylim([GasCoorMin[2], GasCoorMax[2]])
    plt.xlabel('y [ckpc/h]')
    plt.ylabel('z [ckpc/h]')
    # draw_circle = plt.Circle(( CircleCenter[1],  CircleCenter[2]), 150, fill=False)
    # plt.gcf().gca().add_artist(draw_circle)
    plt.colorbar(image)
    plt.gca().set_aspect(1)
    plt.savefig(savepath+'/{}.png'.format(i))

def subhaloID2subfindID (subhaloID,sublink=sublink):
    subhaloID_array=sublink['SubhaloID']
    subfindIndex=np.where(subhaloID_array==subhaloID)
    subfindID=sublink['SubfindID'][subfindIndex]

    return subfindID

def DescendantSubfindID(subfindID,snapshot_num,sublink=sublink):
    tree_fields=['SubfindID','SnapNum','DescendantID']
    subfindid_snapnum_tree=il.sublink.loadTree(basePath, snapshot_num,subfindID, fields=tree_fields, onlyMPB=True)
    if subfindid_snapnum_tree['DescendantID'][0]==-1:return -1
    descendant_subfind_id=subhaloID2subfindID(subfindid_snapnum_tree['DescendantID'][0],sublink)[0]
    return descendant_subfind_id

def batch_distance_calculation(pos1, pos2, boxsize=None):
    # 计算每个点与目标点的距离
    dxs = np.abs(pos1[:, 0] - pos2[0])
    dys = np.abs(pos1[:, 1] - pos2[1])
    dzs = np.abs(pos1[:, 2] - pos2[2])

    # 考虑周期性边界条件
    if boxsize is not None:
        dxs2 = np.abs(boxsize - dxs)
        dys2 = np.abs(boxsize - dys)
        dzs2 = np.abs(boxsize - dzs)

        # 在周期性边界条件下，选择最短的距离
        dxs = np.minimum(dxs, dxs2)
        dys = np.minimum(dys, dys2)
        dzs = np.minimum(dzs, dzs2)

    # 计算三维空间点之间的欧氏距离
    distances = np.linalg.norm(np.vstack([dxs, dys, dzs]).T, axis=1)

    return distances


def Radius90(Particles_Coordinates,Particles_Masses,Center):
    boxsize=51700
    Particles_Distance=batch_distance_calculation(Particles_Coordinates,Center,boxsize)
    Particles_Distance_Sorted=np.sort(Particles_Distance)
    Particles_Masses_Sorted=Particles_Masses[np.argsort(Particles_Distance)]
    Masses_Cumulative=np.cumsum(Particles_Masses_Sorted)
    Masses_Cumulative_Normalized=Masses_Cumulative/Masses_Cumulative[-1]
    Radius90=Particles_Distance_Sorted[np.where(Masses_Cumulative_Normalized>=0.9)[0][0]]
    Radius50=Particles_Distance_Sorted[np.where(Masses_Cumulative_Normalized>=0.5)[0][0]]
    return Radius90,Radius50

def Radius90_Process_And_Save(Subhalo_Index):
    Subhalo_Dict=il.groupcat.loadSingle(basePath, Snapshot_End, -1,Subhalo_Index)
    Subhalo_Center=Subhalo_Dict['SubhaloPos']/h
    Subhalo_Mass=Subhalo_Dict['SubhaloMass']/h
    if Subhalo_Mass>10000:
        HI_Radius90=-1
        Stars_Radius90=-1
        File=dict()
        File['Subhalo_Index']=Subhalo_Index
        File['Subhalo_MassType']=Subhalo_Dict['SubhaloMassType']/h
        File['HI_Radius90']=HI_Radius90
        File['HI_Radius50']=HI_Radius90
        File['Stars_Radius50']=Stars_Radius90
        File['Stars_Radius90']=Stars_Radius90
        hdf5_filename=str(Subhalo_Index)+'.hdf5'
        hdf5_path=os.path.join('./Radius90',hdf5_filename)
        with h5py.File(hdf5_path,'w') as f:
            for key in File.keys():
                f[key]=File[key]
        return 0

    Gas_Dict=il.snapshot.loadSubhalo(basePath, Snapshot_End, Subhalo_Index,'gas',fields=Gas_Fields)
    if len(Gas_Dict)==1:
        HI_Radius90=0
        Stars_Radius90=0
        File=dict()
        File['Subhalo_Index']=Subhalo_Index
        File['Subhalo_MassType']=Subhalo_Dict['SubhaloMassType']/h
        File['HI_Radius90']=HI_Radius90
        File['HI_Radius50']=HI_Radius90
        File['Stars_Radius50']=Stars_Radius90
        File['Stars_Radius90']=Stars_Radius90
        hdf5_filename=str(Subhalo_Index)+'.hdf5'
        hdf5_path=os.path.join('./Radius90',hdf5_filename)
        with h5py.File(hdf5_path,'w') as f:
            for key in File.keys():
                f[key]=File[key]
        return 0
    Gas_Coordinates=Gas_Dict['Coordinates']/h
    Gas_Masses=Gas_Dict['Masses']/h
    Gas_NeutralHydrogenAbundance=Gas_Dict['NeutralHydrogenAbundance']
    Gas_NeutralHydrogen_Masses=Gas_Masses*Gas_NeutralHydrogenAbundance



    HI_Radius90=Radius90(Gas_Coordinates,Gas_NeutralHydrogen_Masses,Subhalo_Center)
    del Gas_Dict,Gas_Coordinates,Gas_Masses,Gas_NeutralHydrogenAbundance,Gas_NeutralHydrogen_Masses

    Stars_Dict=il.snapshot.loadSubhalo(basePath, Snapshot_End, Subhalo_Index,'stars',fields=Stars_Fields)
    Stars_Coordinates=Stars_Dict['Coordinates']/h
    Stars_Masses=Stars_Dict['Masses']/h
    Stars_Radius90=Radius90(Stars_Coordinates,Stars_Masses,Subhalo_Center)
    del Stars_Dict,Stars_Coordinates,Stars_Masses

    File=dict()
    File['Subhalo_Index']=Subhalo_Index
    File['Subhalo_MassType']=Subhalo_Dict['SubhaloMassType']/h
    File['HI_Radius90']=HI_Radius90[0]
    File['HI_Radius50']=HI_Radius90[1]
    File['Stars_Radius90']=Stars_Radius90[0]
    File['Stars_Radius50']=Stars_Radius90[1]


    hdf5_filename=str(Subhalo_Index)+'.hdf5'
    hdf5_path=os.path.join('./Radius90',hdf5_filename)
    with h5py.File(hdf5_path,'w') as f:
        for key in File.keys():
            f[key]=File[key]
    return 0


def process_subhaloindex(i,Subhalo_Indices_Chunks):
    Subhalo_Indices=Subhalo_Indices_Chunks[i]
    for Subhalo_Index in tqdm(Subhalo_Indices):
        Radius90_Process_And_Save(Subhalo_Index)

def chunks(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def mass_within_radius(Radius, Masses, Coordinates, CenterOfMass):
    # 计算每个恒星到质心的距离
    dist_to_com = np.sqrt(np.sum((Coordinates - CenterOfMass)**2, axis=1))

    # 找到距离质心小于R_200的所有恒星的下标
    within_radius_idx = np.where(dist_to_com <= Radius)[0]

    # 计算在R_200范围内的恒星总质量
    total_mass_within_radius = np.sum(Masses[within_radius_idx])

    return total_mass_within_radius

def baryon_mass_in_radius(Center, Radius, Coordinates, Masses):
    dist_to_com = np.sqrt(np.sum((Coordinates - Center)**2, axis=1))
    within_radius_idx = np.where(dist_to_com <= Radius)[0]
    total_mass_within_radius = np.sum(Masses[within_radius_idx])
    return total_mass_within_radius

def drawsubhalo(i,sp,savepath,basePath):
    index=np.where(sorted_index_array==i)

    gas_data=il.snapshot.loadSubhalo(basePath, sp, i, 'gas', fields=['Coordinates','NeutralHydrogenAbundance','Masses'])
    if len(gas_data)==1:return -1

    stars_mass=il.snapshot.loadSubhalo(basePath,current_snap_num,i,'stars',fields=['Masses'])

    cold_gas_mass =Cold_Gas_Mass(sp,i)

    GasCoordinates = gas_data['Coordinates']
    GasCoordinates=np.mod(GasCoordinates,35000)
    GasAbundance = gas_data['NeutralHydrogenAbundance']
    Masses = gas_data['Masses']
    NeutralHydrogenAbundance=Masses*GasAbundance

    GasCoorMin = GasCoordinates.min(axis=0)
    GasCoorMax = GasCoordinates.max(axis=0)
    bars=2*(GasCoorMax-GasCoorMin)
    if (bars[0]>8000)or (bars[2]>8000):return -1
    plt.figure(figsize=(20, 20), dpi=400)

    plt.style.use("dark_background")

    h, _, _, image = plt.hist2d(GasCoordinates[:, 0],
                            GasCoordinates[:, 2],
                            weights=NeutralHydrogenAbundance,

                            norm=mpl.colors.LogNorm(),
                            bins=(bars[0],bars[2]))
    plt.xlim([GasCoorMin[0],GasCoorMax[0]])
    plt.ylim([GasCoorMin[2],GasCoorMax[2]])
    plt.xlabel('x [ckpc/h]')
    plt.ylabel('z [ckpc/h]')
    plt.text(0.65*(GasCoorMax[0]-GasCoorMin[0])+GasCoorMin[0],0.82*(GasCoorMax[2]-GasCoorMin[2])+GasCoorMin[2],'HI_R_90={:.2f}kpc\nM_Stars={:.2f}E10Ms\nM_Galaxy={:.2f}E10Ms\nM_Gas={:.2f}E10Ms\nM_CG={:.2f}E10Ms\nM_HI={:.2f}E10Ms'.format(sorted_R_90_array[np.where(sorted_index_array==i)[0]][0],stars_mass.sum(),subhalo_mass[i],Masses.sum(),cold_gas_mass,NeutralHydrogenAbundance.sum()),color='white',fontsize=22)

    draw_circle = plt.Circle((subhalo_pos[i][0], subhalo_pos[i][2]),sorted_R_90_array[index][0] , fill=False)
    plt.gcf().gca().add_artist(draw_circle)
    plt.colorbar(image)
    plt.gca().set_aspect(1)
    plt.savefig(savepath+'/{}.png'.format(i))
    plt.close('all')
    plt.clf()
    del index, GasCoorMin,GasCoorMax,gas_data,stars_mass,cold_gas_mass,GasCoordinates,GasAbundance,Masses,NeutralHydrogenAbundance,bars,image,draw_circle
    gc.collect()
    return 0


def calculate_R90(Rmin, Rmax, MassRatio,Masses,Coordinates,CenterOfMass):
    # 定义精度和最大迭代次数
    precision = 1
    max_iter = 1000

    # 初始化左右边界
    left, right = Rmin, Rmax

    # 迭代计算R_90
    for i in range(max_iter):
        # 计算中间点
        mid = (left + right) / 2

        # 计算该点的质量占比
        ratio =  mass_ratio_within_radius(mid, Masses, Coordinates, CenterOfMass)

        # 如果质量占比小于90%，则将左边界更新为mid
        if ratio < 0.9:
            left = mid
        # 否则将右边界更新为mid
        else:
            right = mid

        # 如果当前的区间长度小于给定的精度，或者达到最大迭代次数，就停止迭代
        if abs(right - left) < precision or i == max_iter - 1:
            R90 = (left + right) / 2
            break

    return R90

def mass_ratio_within_radius(Radius, Masses, Coordinates, CenterOfMass):
    # 计算每个恒星到质心的距离
    dist_to_com = np.sqrt(np.sum((Coordinates - CenterOfMass)**2, axis=1))

    # 找到距离质心小于R_200的所有恒星的下标
    within_radius_idx = np.where(dist_to_com <= Radius)[0]

    # 计算在R_200范围内的恒星总质量
    total_mass_within_radius = np.sum(Masses[within_radius_idx])

    # 计算星系总质量
    total_mass = np.sum(Masses)

    # 计算质量比值
    mass_ratio = total_mass_within_radius / total_mass

    return mass_ratio

def Binary_DataIter(Data_Set1,Data_Set2):
    for i in range(0,min(len(Data_Set1),len(Data_Set2))):
        yield Data_Set1[i],Data_Set2[i]

def R_200_Calculation(snapshot_num,subhalo_id,Critial_Density):

    subhalos_data=il.groupcat.loadSubhalos(basePath,snapshot_num,fields=subhalo_fields)
    subhalo_position=subhalos_data['SubhaloPos'][subhalo_id]
    subhalo_velocity=subhalos_data['SubhaloVel'][subhalo_id]
    subhalo_spin=subhalos_data['SubhaloSpin'][subhalo_id]
    subhalo_mass=subhalos_data['SubhaloMass'][subhalo_id]


    gas_data_tensor=il.snapshot.loadSubhalo(basePath,snapshot_num,subhalo_id,0,fields=gas_fields)
    dm_data_tensor=il.snapshot.loadSubhalo(basePath,snapshot_num,subhalo_id,1,fields=dm_fields)
    stars_data_tensor=il.snapshot.loadSubhalo(basePath,snapshot_num,subhalo_id,4,fields=stars_fields)
    bh_data_tensor=il.snapshot.loadSubhalo(basePath,snapshot_num,subhalo_id,5,fields=bh_fields)


    precision = 1
    max_iter = 100
    left, right = 0, 1000

    for i in range(max_iter):

        mid = (left + right) / 2

        mass_in=total_mass_in_radius(gas_data_tensor,dm_data_tensor,stars_data_tensor,bh_data_tensor,subhalo_position,mid)

        density=3*mass_in/(4*3.14*mid**3)

        if density > 200*Critial_Density:
            left = mid

        else:
            right = mid


        if abs(right - left) < precision or i == max_iter - 1:
            R_200 = (left + right) / 2
            break

    return R_200,subhalo_mass

def baryon_mass_in_radius(Center,Radius,Coordinates,Masses):
    dist_to_com = np.sqrt(np.sum((Coordinates - Center)**2, axis=1))
    within_radius_idx = np.where(dist_to_com <= Radius)[0]
    total_mass_within_radius = np.sum(Masses[within_radius_idx])
    return total_mass_within_radius

def dm_mass_in_radius(Center,Radius,Coordinates):
    dist_to_com = np.sqrt(np.sum((Coordinates - Center)**2, axis=1))
    within_radius_idx = np.where(dist_to_com <= Radius)[0]
    total_mass_within_radius = m_dm*within_radius_idx.shape[0]*h
    return total_mass_within_radius

def total_mass_in_radius(gas_data_tensor,dm_data_tensor,stars_data_tensor,bh_data_tensor,Center,Radius):

    gas_mass_in=baryon_mass_in_radius(Center,Radius,gas_data_tensor['Coordinates'],gas_data_tensor['Masses'])
    dm_mass_in=dm_mass_in_radius(Center,Radius,dm_data_tensor['Coordinates'])
    stars_mass_in=baryon_mass_in_radius(Center,Radius,stars_data_tensor['Coordinates'],stars_data_tensor['Masses'])
    bh_mass_in=baryon_mass_in_radius(Center,Radius,bh_data_tensor['Coordinates'],bh_data_tensor['Masses'])

    total_mass_in=gas_mass_in+dm_mass_in+stars_mass_in+bh_mass_in

    return total_mass_in

def cal_dm_angular_momentum(Position,Velocity,DM_Data_Tensor,R_virial):


    DM_Coordinates,DM_Velocities=DM_Data_Tensor['Coordinates'],DM_Data_Tensor['Velocities']

    dm_relative_velocities_tensor=DM_Velocities-Velocity
    dm_relative_coordinates_tensor=DM_Coordinates-Position
    dm_particle_number=dm_relative_coordinates_tensor.shape[0]
    dm_particle_distance=calculate_distance(dm_relative_coordinates_tensor)

    cross_product=np.cross(dm_relative_coordinates_tensor,dm_relative_velocities_tensor)
    index=np.where(dm_particle_distance>R_virial)[0]

    for i in index:
        cross_product[i]=0


    dm_angular_momentum=(cross_product/dm_particle_number).sum(axis=0)

    return dm_angular_momentum

def calculate_distance(array):
    squared=np.square(array)

    sum=squared.sum(axis=1)
    sum=np.sqrt(sum)
    return sum

def cal_baryon_angular_momentum(Position,Velocity,Baryon_Data_Tensor,R_virial):

    Baryon_Coordinates,Baryon_Velocities,Baryon_Masses=Baryon_Data_Tensor['Coordinates'],Baryon_Data_Tensor['Velocities'],Baryon_Data_Tensor['Masses']

    baryon_relative_velocities_tensor=Baryon_Velocities-Velocity
    baryon_relative_coordinates_tensor=Baryon_Coordinates-Position
    baryon_particle_distance=calculate_distance(baryon_relative_coordinates_tensor)

    cross_product=np.cross(baryon_relative_coordinates_tensor,baryon_relative_velocities_tensor)

    index=np.where(baryon_particle_distance>R_virial)[0]

    for i in index:
        cross_product[i]=0

    angular_momentum=(Baryon_Masses.reshape(-1,1)*cross_product).sum(axis=0)
    weight_angular_momentum=(angular_momentum/Baryon_Masses.sum())

    return weight_angular_momentum


def HI_Radius_90_Calculation_At_99(Subhalo_Id):
    max_iter=100
    precension=1
    left = 0
    right=500

    subhalo_data_dict = il.groupcat.loadSingle(basePath,99,-1,Subhalo_Id)

    subhalo_comoving_position_floatarray = subhalo_data_dict['SubhaloPos']

    gas_data_dict=il.snapshot.loadSubhalo(basePath,99,Subhalo_Id,'gas',fields=['Coordinates','Masses','NeutralHydrogenAbundance'])

    if len(gas_data_dict)==1: return 0

    gas_data_dict['Masses']=gas_data_dict['Masses']*gas_data_dict['NeutralHydrogenAbundance']

    gas_mass_90=gas_data_dict['Masses'].sum()*0.9

    for epoch in range(max_iter):
        mid=(left+right)/2

        mass_in=Baryon_Mass_Within_Given_Radius(subhalo_comoving_position_floatarray,mid,gas_data_dict['Coordinates'],gas_data_dict['Masses'])

        if mass_in > gas_mass_90: right=mid

        else: left=mid

        if abs(right-left)<precension or epoch == max_iter:
            R_90_Com_Float = (left + right) / 2
            break

    del max_iter,precension,left,right,subhalo_data_dict,subhalo_comoving_position_floatarray,gas_data_dict,gas_mass_90

    return R_90_Com_Float

def Stellar_Radius_90_Calculation_At_99(Subhalo_Id):
    max_iter=100
    precension=1
    left = 0
    right=500

    subhalo_data_dict = il.groupcat.loadSingle(basePath,99,-1,Subhalo_Id)

    subhalo_comoving_position_floatarray = subhalo_data_dict['SubhaloPos']

    stellar_data_dict=il.snapshot.loadSubhalo(basePath,99,Subhalo_Id,'stellar',fields=['Coordinates','Masses'])

    if len(stellar_data_dict)==1: return 0

    stellar_mass_90=stellar_data_dict['Masses'].sum()*0.9

    for epoch in range(max_iter):

        mid=(left+right)/2
        mass_in=Baryon_Mass_Within_Given_Radius(subhalo_comoving_position_floatarray,mid,stellar_data_dict['Coordinates'],stellar_data_dict['Masses'])

        if mass_in > stellar_mass_90: right=mid

        else: left=mid

        if abs(right-left)<precension or epoch == max_iter:
            R_90_Com_Float = (left + right) / 2
            break

    del max_iter,precension,left,right,subhalo_data_dict,subhalo_comoving_position_floatarray,stellar_data_dict,stellar_mass_90

    return R_90_Com_Float



def Baryon_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,Particle_Coordinates_FloatArray,Particle_Masses_FloatArray):
    dist_to_com_floatarray = np.sqrt(np.sum((Particle_Coordinates_FloatArray - Galaxy_Center_FloatArray)**2, axis=1))
    within_radius_indices_intarray = np.where(dist_to_com_floatarray <= Mass_Radius_Float)[0]
    Total_Mass_Within_Radius_Float = np.sum(Particle_Masses_FloatArray[within_radius_indices_intarray])

    del dist_to_com_floatarray,within_radius_indices_intarray

    return Total_Mass_Within_Radius_Float

def Total_Mass_Type_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,Data_Dict):

    mass_within_radius=Baryon_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,Data_Dict['Coordinates'],Data_Dict['Masses'])

    return mass_within_radius

def Calculate_Radius_Ratio(Subhalo_Id):
    R_HI=HI_Radius_90_Calculation_At_99(Subhalo_Id)
    R_Stellar=Stellar_Radius_90_Calculation_At_99(Subhalo_Id)

    if R_HI * R_Stellar:
        return R_HI/R_Stellar

    else: return -1

def R_200_Calculation(Snapshot_Num, Subhalo_Id):
    subhalo_data_dict = il.groupcat.loadSingle(basePath, Snapshot_Num, -1, Subhalo_Id)

    subhalo_comoving_position_floatarray = subhalo_data_dict["SubhaloPos"]

    gas_data_dict = il.snapshot.loadSubhalo(
        basePath, Snapshot_Num, Subhalo_Id, "gas", fields=["Coordinates", "Masses"]
    )
    darkmatter_data_dict = il.snapshot.loadSubhalo(
        basePath, Snapshot_Num, Subhalo_Id, "darkmatter", fields=["Coordinates"]
    )
    stellar_data_dict = il.snapshot.loadSubhalo(
        basePath, Snapshot_Num, Subhalo_Id, "stellar", fields=["Coordinates", "Masses"]
    )
    blackholes_data_dict = il.snapshot.loadSubhalo(
        basePath,
        Snapshot_Num,
        Subhalo_Id,
        "blackholes",
        fields=["Coordinates", "Masses"],
    )

    precision = 1
    max_iter = 100
    left, right = 0, 500

    header_dict = il.groupcat.loadHeader(basePath, Snapshot_Num)

    scale_factor = header_dict["Time"]

    for i in range(max_iter):
        mid = (left + right) / 2

        mass_in_radius = Total_Mass_Within_Given_Radius(
            subhalo_comoving_position_floatarray,
            mid,
            gas_data_dict,
            darkmatter_data_dict,
            stellar_data_dict,
            blackholes_data_dict,
        )

        mass_200 = Cal_M_From_R_Phy(mid * scale_factor, header_dict)

        if mass_in_radius > mass_200:
            left = mid

        else:
            right = mid

        if abs(right - left) < precision or i == max_iter - 1:
            R_200_Com_Float = (left + right) / 2
            break

    R_200_Phy_Float = R_200_Com_Float * scale_factor

    M_200_Float = Cal_M_From_R_Phy(R_200_Phy_Float, header_dict)

    del (
        subhalo_data_dict,
        subhalo_comoving_position_floatarray,
        gas_data_dict,
        darkmatter_data_dict,
        stellar_data_dict,
        blackholes_data_dict,
        precision,
        max_iter,
        left,
        right,
        header_dict,
        scale_factor,
        mid,
        mass_in_radius,
        mass_200,
    )

    gc.collect()
    return R_200_Phy_Float, M_200_Float


def Cal_M_From_R_Phy(R_200_Float, Header_Dict):
    Omega_m = Header_Dict["Omega0"]
    Omega_Lambda = Header_Dict["OmegaLambda"]
    z = Header_Dict["Redshift"]
    scale_factor = Header_Dict["Time"]

    E_z = Omega_m * (1 + z) ** 3 + Omega_Lambda
    M_200_Float = 200 * Critial_Density_h * E_z * (R_200_Float**3) * 4 / 3 * 3.14

    del Omega_m, Omega_Lambda, z, scale_factor, E_z

    return M_200_Float


def Baryon_Mass_Within_Given_Radius(
    Galaxy_Center_FloatArray,
    Mass_Radius_Float,
    Particle_Coordinates_FloatArray,
    Particle_Masses_FloatArray,
):
    dist_to_com_floatarray = np.sqrt(
        np.sum(
            (Particle_Coordinates_FloatArray - Galaxy_Center_FloatArray) ** 2, axis=1
        )
    )
    within_radius_indices_intarray = np.where(
        dist_to_com_floatarray <= Mass_Radius_Float
    )[0]
    Total_Mass_Within_Radius_Float = np.sum(
        Particle_Masses_FloatArray[within_radius_indices_intarray]
    )

    del dist_to_com_floatarray, within_radius_indices_intarray

    return Total_Mass_Within_Radius_Float


def DarkMatter_Mass_Within_Given_Radius(
    Galaxy_Center_FloatArray, Mass_Radius_Float, Particle_Coordinates_FloatArray
):
    dist_to_com_floatarray = np.sqrt(
        np.sum(
            (Particle_Coordinates_FloatArray - Galaxy_Center_FloatArray) ** 2, axis=1
        )
    )
    within_radius_indices_intarray = np.where(
        dist_to_com_floatarray <= Mass_Radius_Float
    )[0]
    Total_Mass_Within_Radius_Float = m_dm_h * within_radius_indices_intarray.shape[0]

    del dist_to_com_floatarray, within_radius_indices_intarray

    return Total_Mass_Within_Radius_Float


def Total_Mass_Within_Given_Radius(
    Galaxy_Center_FloatArray,
    Mass_Radius_Float,
    Gas_Data_Dict,
    DarkMatter_Data_Dict,
    Stellar_Data_Dict,
    BlackHole_Data_Dict,
):
    darkmatter_mass_within_radius = DarkMatter_Mass_Within_Given_Radius(
        Galaxy_Center_FloatArray, Mass_Radius_Float, DarkMatter_Data_Dict
    )

    if len(Gas_Data_Dict) == 1:
        gas_mass_within_radius = 0
    else:
        gas_mass_within_radius = Baryon_Mass_Within_Given_Radius(
            Galaxy_Center_FloatArray,
            Mass_Radius_Float,
            Gas_Data_Dict["Coordinates"],
            Gas_Data_Dict["Masses"],
        )

    if len(Stellar_Data_Dict) == 1:
        stellar_mass_within_radius = 0
    else:
        stellar_mass_within_radius = Baryon_Mass_Within_Given_Radius(
            Galaxy_Center_FloatArray,
            Mass_Radius_Float,
            Stellar_Data_Dict["Coordinates"],
            Stellar_Data_Dict["Masses"],
        )

    if len(BlackHole_Data_Dict) == 1:
        blackholes_mass_within_radius = 0
    else:
        blackholes_mass_within_radius = Baryon_Mass_Within_Given_Radius(
            Galaxy_Center_FloatArray,
            Mass_Radius_Float,
            BlackHole_Data_Dict["Coordinates"],
            BlackHole_Data_Dict["Masses"],
        )

    Total_Mass_Within = (
        gas_mass_within_radius
        + darkmatter_mass_within_radius
        + stellar_mass_within_radius
        + blackholes_mass_within_radius
    )

    del (
        gas_mass_within_radius,
        darkmatter_mass_within_radius,
        stellar_mass_within_radius,
        blackholes_mass_within_radius,
    )

    return Total_Mass_Within


def Cal_Subhalo_Id_at_Z_from_99(Subhalo_Id_0, Z):  # Z is the snapshot
    subfindid_snapnum_tree = il.sublink.loadTree(
        basePath,
        99,
        Subhalo_Id_0,
        fields=["SubfindID", "SnapNum", "Mass"],
        onlyMPB=True,
    )

    if subfindid_snapnum_tree == None:
        return -1

    subfind_id = subfindid_snapnum_tree["SubfindID"]
    snapshot_num = subfindid_snapnum_tree["SnapNum"]
    index = np.where(snapshot_num == Z)[0]
    if index.shape[0] == 0:
        return -1
    return subfind_id[index][0]

def Cal_M_From_R_Phy(R_200_Float,Header_Dict):

    Omega_m=Header_Dict['Omega0']
    Omega_Lambda=Header_Dict['OmegaLambda']
    z=Header_Dict['Redshift']
    scale_factor=Header_Dict['Time']

    E_z=Omega_m*(1+z)**3+Omega_Lambda
    M_200_Float=200*Critial_Density_h*E_z*(R_200_Float**3)*4/3*3.14

    del Omega_m,Omega_Lambda,z,scale_factor,E_z

    return M_200_Float

def Baryon_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,Particle_Coordinates_FloatArray,Particle_Masses_FloatArray):
    dist_to_com_floatarray = np.sqrt(np.sum((Particle_Coordinates_FloatArray - Galaxy_Center_FloatArray)**2, axis=1))
    within_radius_indices_intarray = np.where(dist_to_com_floatarray <= Mass_Radius_Float)[0]
    Total_Mass_Within_Radius_Float = np.sum(Particle_Masses_FloatArray[within_radius_indices_intarray])

    del dist_to_com_floatarray,within_radius_indices_intarray

    return Total_Mass_Within_Radius_Float

def DarkMatter_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,Particle_Coordinates_FloatArray):
    dist_to_com_floatarray = np.sqrt(np.sum((Particle_Coordinates_FloatArray - Galaxy_Center_FloatArray)**2, axis=1))
    within_radius_indices_intarray = np.where(dist_to_com_floatarray <= Mass_Radius_Float)[0]
    Total_Mass_Within_Radius_Float = m_dm_h*within_radius_indices_intarray.shape[0]

    del dist_to_com_floatarray,within_radius_indices_intarray

    return Total_Mass_Within_Radius_Float

def Total_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,Gas_Data_Dict,DarkMatter_Data_Dict,Stellar_Data_Dict,BlackHole_Data_Dict):

    darkmatter_mass_within_radius=DarkMatter_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,DarkMatter_Data_Dict)

    if len(Gas_Data_Dict)==1:gas_mass_within_radius=0
    else: gas_mass_within_radius=Baryon_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,Gas_Data_Dict['Coordinates'],Gas_Data_Dict['Masses'])

    if len(Stellar_Data_Dict)==1:stellar_mass_within_radius=0
    else: stellar_mass_within_radius=Baryon_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,Stellar_Data_Dict['Coordinates'],Stellar_Data_Dict['Masses'])


    if len(BlackHole_Data_Dict)==1:blackholes_mass_within_radius=0
    else: blackholes_mass_within_radius=Baryon_Mass_Within_Given_Radius(Galaxy_Center_FloatArray,Mass_Radius_Float,BlackHole_Data_Dict['Coordinates'],BlackHole_Data_Dict['Masses'])

    Total_Mass_Within=gas_mass_within_radius+darkmatter_mass_within_radius+stellar_mass_within_radius+blackholes_mass_within_radius

    del gas_mass_within_radius,darkmatter_mass_within_radius,stellar_mass_within_radius,blackholes_mass_within_radius

    return Total_Mass_Within

def Cal_Subhalo_Id_at_Z_from_99(Subhalo_Id_0,Z):#Z is the snapshot
    subfindid_snapnum_tree=il.sublink.loadTree(basePath, 99, Subhalo_Id_0, fields=['SubfindID','SnapNum','Mass'], onlyMPB=True)

    if subfindid_snapnum_tree == None:return -1

    subfind_id=subfindid_snapnum_tree['SubfindID']
    snapshot_num=subfindid_snapnum_tree['SnapNum']
    index=np.where(snapshot_num==Z)[0]
    if index.shape[0] ==0:return -1
    return subfind_id[index][0]

def Velocity_Virial_Calculation(R_200,M_200):
    Unit=43007.1
    V=np.sqrt(M_200/R_200*Unit)
    return V

def Cal_Total_Spin(spin_array):
    len=0.0
    for axis_value in spin_array:
        len=len+axis_value**2

    return np.sqrt(len)

def R_200_Calculation(Snapshot_Num,Subhalo_Id):

    subhalo_data_dict=il.groupcat.loadSingle(basePath,Snapshot_Num,-1,Subhalo_Id)

    subhalo_comoving_position_floatarray=subhalo_data_dict['SubhaloPos']

    gas_data_dict=il.snapshot.loadSubhalo(basePath,Snapshot_Num,Subhalo_Id,'gas',fields=['Coordinates','Masses'])
    darkmatter_data_dict=il.snapshot.loadSubhalo(basePath,Snapshot_Num,Subhalo_Id,'darkmatter',fields=['Coordinates'])
    stellar_data_dict=il.snapshot.loadSubhalo(basePath,Snapshot_Num,Subhalo_Id,'stellar',fields=['Coordinates','Masses'])
    blackholes_data_dict=il.snapshot.loadSubhalo(basePath,Snapshot_Num,Subhalo_Id,'blackholes',fields=['Coordinates','Masses'])

    precision = 1
    max_iter = 100
    left, right = 0, 500

    header_dict=il.groupcat.loadHeader(basePath,Snapshot_Num)

    scale_factor=header_dict['Time']

    for i in range(max_iter):

        mid = (left + right) / 2

        mass_in_radius=Total_Mass_Within_Given_Radius(subhalo_comoving_position_floatarray,mid,gas_data_dict,darkmatter_data_dict,stellar_data_dict,blackholes_data_dict)

        mass_200=Cal_M_From_R_Phy(mid*scale_factor,header_dict)

        if mass_in_radius > mass_200:
            left = mid

        else:
            right = mid


        if abs(right - left) < precision or i == max_iter - 1:
            R_200_Com_Float = (left + right) / 2
            break

    R_200_Phy_Float=R_200_Com_Float*scale_factor

    M_200_Float=Cal_M_From_R_Phy(R_200_Phy_Float,header_dict)

    del subhalo_data_dict,subhalo_comoving_position_floatarray, gas_data_dict,darkmatter_data_dict,stellar_data_dict,blackholes_data_dict,precision,max_iter, left, right,header_dict,scale_factor,mid,mass_in_radius,mass_200

    gc.collect()
    return R_200_Phy_Float,M_200_Float

def Spin_Parameter_Calculation(Snapshot_Num,Subhalo_Id):


    header_dict=il.groupcat.loadHeader(basePath,Snapshot_Num)
    scale_factor=header_dict['Time']

    R_200, M_200 = R_200_Calculation(Snapshot_Num, Subhalo_Id)
    V_200 = Velocity_Virial_Calculation(R_200, M_200)

    subhalo_data_dict=il.groupcat.loadSingle(basePath,Snapshot_Num,-1,Subhalo_Id)

    total_spin=Cal_Total_Spin(subhalo_data_dict['SubhaloSpin'])

    Spin_Parameter=total_spin*scale_factor/(1.414*V_200*R_200)

    del R_200,M_200,header_dict,scale_factor,V_200,subhalo_data_dict,total_spin

    gc.collect()
    return Spin_Parameter

def Cal_M_200_From_R_200(R_200,Comoving_Critial_Density):

    M_200=200*Comoving_Critial_Density*4/3*3.14*R_200**3
    return M_200

def R_200_Calculation(snapshot_num,subhalos_data,subhalo_id,Comoving_Critial_Density):

    subhalo_comoving_position=subhalos_data['SubhaloPos'][subhalo_id]
    subhalo_comoving_velocity=subhalos_data['SubhaloVel'][subhalo_id]
    subhalo_spin=subhalos_data['SubhaloSpin'][subhalo_id]
    subhalo_mass=subhalos_data['SubhaloMass'][subhalo_id]

    gas_data_tensor=il.snapshot.loadSubhalo(basePath,snapshot_num,subhalo_id,0,fields=gas_fields)
    dm_data_tensor=il.snapshot.loadSubhalo(basePath,snapshot_num,subhalo_id,1,fields=dm_fields)
    stars_data_tensor=il.snapshot.loadSubhalo(basePath,snapshot_num,subhalo_id,4,fields=stars_fields)
    bh_data_tensor=il.snapshot.loadSubhalo(basePath,snapshot_num,subhalo_id,5,fields=bh_fields)


    precision = 0.1
    max_iter = 100
    left, right = 0, 1000

    for i in range(max_iter):

        mid = (left + right) / 2

        mass_in=total_mass_in_radius(gas_data_tensor,dm_data_tensor,stars_data_tensor,bh_data_tensor,subhalo_comoving_position,mid)

        density=3*mass_in/(4*3.14*mid**3)

        if density > 200*Comoving_Critial_Density:
            left = mid

        else:
            right = mid


        if abs(right - left) < precision or i == max_iter - 1:
            R_200 = (left + right) / 2
            break

    M_200=Cal_M_200_From_R_200(R_200,Comoving_Critial_Density)

    return R_200,M_200

def baryon_mass_in_radius(Center,Radius,Coordinates,Masses):
    dist_to_com = np.sqrt(np.sum((Coordinates - Center)**2, axis=1))
    within_radius_idx = np.where(dist_to_com <= Radius)[0]
    total_mass_within_radius = np.sum(Masses[within_radius_idx])
    return total_mass_within_radius

def dm_mass_in_radius(Center,Radius,Coordinates):
    dist_to_com = np.sqrt(np.sum((Coordinates - Center)**2, axis=1))
    within_radius_idx = np.where(dist_to_com <= Radius)[0]
    total_mass_within_radius = m_dm*within_radius_idx.shape[0]*h
    return total_mass_within_radius

def total_mass_in_radius(gas_data_tensor,dm_data_tensor,stars_data_tensor,bh_data_tensor,Center,Radius):

    gas_mass_in=baryon_mass_in_radius(Center,Radius,gas_data_tensor['Coordinates'],gas_data_tensor['Masses'])
    dm_mass_in=dm_mass_in_radius(Center,Radius,dm_data_tensor['Coordinates'])
    stars_mass_in=baryon_mass_in_radius(Center,Radius,stars_data_tensor['Coordinates'],stars_data_tensor['Masses'])
    bh_mass_in=baryon_mass_in_radius(Center,Radius,bh_data_tensor['Coordinates'],bh_data_tensor['Masses'])

    total_mass_in=gas_mass_in+dm_mass_in+stars_mass_in+bh_mass_in

    return total_mass_in

def Binary_DataIter(Data_Set1,Data_Set2):
    for i in range(0,min(len(Data_Set1),len(Data_Set2))):
        yield Data_Set1[i],Data_Set2[i]

def Velocity_Virial_Calculation(R_200,M_200):
    Unit=4.3E4
    V=np.sqrt(R_200/M_200*Unit)
    return V

def cal_id_at_Z_from_99(Subhalo_Id_0,Z):#Z is the snapshot
    subfindid_snapnum_tree=il.sublink.loadTree(basePath, 99, Subhalo_Id_0, fields=tree_fields, onlyMPB=True)
    subfind_id=subfindid_snapnum_tree['SubfindID']
    snapshot_num=subfindid_snapnum_tree['SnapNum']
    index=np.where(snapshot_num==Z)[0]
    return subfind_id[index][0]

def Cal_Total_Spin(spin_array):
    len=0.0
    for axis_value in spin_array:
        len=len+axis_value**2

    return np.sqrt(len)

def batch_distance_calculation(pos1, pos2, boxsize=None):
    # 计算每个点与目标点的距离
    dxs = np.abs(pos1[:, 0] - pos2[0])
    dys = np.abs(pos1[:, 1] - pos2[1])
    dzs = np.abs(pos1[:, 2] - pos2[2])

    # 考虑周期性边界条件
    if boxsize is not None:
        dxs2 = np.abs(boxsize - dxs)
        dys2 = np.abs(boxsize - dys)
        dzs2 = np.abs(boxsize - dzs)

        # 在周期性边界条件下，选择最短的距离
        dxs = np.minimum(dxs, dxs2)
        dys = np.minimum(dys, dys2)
        dzs = np.minimum(dzs, dzs2)

    # 计算三维空间点之间的欧氏距离
    distances = np.linalg.norm(np.vstack([dxs, dys, dzs]).T, axis=1)

    return distances


def Radius90(Particles_Coordinates,Particles_Masses,Center):
    boxsize=51700
    Particles_Distance=batch_distance_calculation(Particles_Coordinates,Center,boxsize)
    Particles_Distance_Sorted=np.sort(Particles_Distance)
    Particles_Masses_Sorted=Particles_Masses[np.argsort(Particles_Distance)]
    Masses_Cumulative=np.cumsum(Particles_Masses_Sorted)
    Masses_Cumulative_Normalized=Masses_Cumulative/Masses_Cumulative[-1]
    Radius90=Particles_Distance_Sorted[np.where(Masses_Cumulative_Normalized>=0.9)[0][0]]
    Radius60=Particles_Distance_Sorted[np.where(Masses_Cumulative_Normalized>=0.6)[0][0]]
    Radius50=Particles_Distance_Sorted[np.where(Masses_Cumulative_Normalized>=0.5)[0][0]]
    return Radius90,Radius60,Radius50

def Density_Radius(Particles_Coordinates,Particles_Masses,Subhalo_Center,Ratio):
    boxsize=51700
    Particles_Distance=batch_distance_calculation(Particles_Coordinates,Subhalo_Center,boxsize)
    Particles_Distance_Sorted=np.sort(Particles_Distance)
    Particles_Masses_Sorted=Particles_Masses[np.argsort(Particles_Distance)]
    Masses_Cumulative=np.cumsum(Particles_Masses_Sorted)
    Density_Cumulative=Masses_Cumulative/(4/3*np.pi*np.power(Particles_Distance_Sorted,3))
    Density_Cumulative_Normalized=Density_Cumulative/Critial_Density/200
    Radius=Particles_Distance_Sorted[np.argmax(Density_Cumulative_Normalized):][np.where(Density_Cumulative_Normalized[np.argmax(Density_Cumulative_Normalized):]<Ratio)[0][0]]

    return Radius

def Spin_Calculator(Snapshot,Host_Dict,Satellite_Index):
    Header=il.groupcat.loadHeader(basePath,Snapshot)

    Satellite_Gas_Dict=il.snapshot.loadSubhalo(basePath,Snapshot,Satellite_Index,'gas',fields=['Coordinates','Velocities','Masses'])
    Satellite_Stars_Dict=il.snapshot.loadSubhalo(basePath,Snapshot,Satellite_Index,'stars',fields=['Coordinates','Velocities','Masses'])
    Satellite_BH_Dict=il.snapshot.loadSubhalo(basePath,Snapshot,Satellite_Index,'5',fields=['Coordinates','Velocities','Masses'])
    Satellite_DM_Dict=il.snapshot.loadSubhalo(basePath,Snapshot,Satellite_Index,'dm',fields=['Coordinates','Velocities'])

    Total_Spin=Baryon_Spin_Calculator_From_Dict(Snapshot,Host_Dict,Satellite_Gas_Dict)+Baryon_Spin_Calculator_From_Dict(Snapshot,Host_Dict,Satellite_Stars_Dict)+Baryon_Spin_Calculator_From_Dict(Snapshot,Host_Dict,Satellite_BH_Dict)+DarkMatter_Spin_Calculator(Snapshot,Host_Dict,Satellite_DM_Dict)
    return Total_Spin

def Baryon_Spin_Calculator_From_Dict(Snapshot,Host_Dict,Particle_Dict):
    Header=il.groupcat.loadHeader(basePath,Snapshot)
    scale_factor=Header['Time']

    Host_Center=Host_Dict['SubhaloPos']*scale_factor
    Host_Velocity=Host_Dict['SubhaloVel']

    Particle_Coordinates=Particle_Dict['Coordinates']*scale_factor
    Particle_Velocities=Particle_Dict['Velocities']*np.sqrt(scale_factor)

    Particle_Coordinates=Particle_Coordinates-Host_Center
    Particle_Velocities=Particle_Velocities-Host_Velocity

    Particle_Angular_Momentum=np.cross(Particle_Coordinates,Particle_Velocities)
    Particle_Angular_Momentum=Particle_Angular_Momentum*Particle_Dict['Masses'].reshape(-1,1)
    Particle_Angular_Momentum_Sum=np.sum(Particle_Angular_Momentum,axis=0)

    return Particle_Angular_Momentum_Sum

def DarkMatter_Spin_Calculator(Snapshot,Host_Dict,Particle_Dict):
    Header=il.groupcat.loadHeader(basePath,Snapshot)
    scale_factor=Header['Time']

    Host_Center=Host_Dict['SubhaloPos']*scale_factor
    Host_Velocity=Host_Dict['SubhaloVel']

    Particle_Coordinates=Particle_Dict['Coordinates']*scale_factor
    Particle_Velocities=Particle_Dict['Velocities']*np.sqrt(scale_factor)

    Particle_Coordinates=Particle_Coordinates-Host_Center
    Particle_Velocities=Particle_Velocities-Host_Velocity

    Particle_Angular_Momentum=np.cross(Particle_Coordinates,Particle_Velocities)
    Particle_Angular_Momentum=Particle_Angular_Momentum*m_dm_h
    Particle_Angular_Momentum_Sum=np.sum(Particle_Angular_Momentum,axis=0)

    return Particle_Angular_Momentum_Sum

def Cold_Gas_Mass_Calculator(Snapshot, Subhalo_Index):
    Gas_Dict=il.snapshot.loadSubhalo(basePath,Snapshot,Subhalo_Index,'gas',fields=['Masses','InternalEnergy','ElectronAbundance'])
    if len(Gas_Dict)==1:return 0,0
    Gas_Masses=Gas_Dict['Masses'].astype(np.float64)
    Gas_Internal_Energy=Gas_Dict['InternalEnergy'].astype(np.float64)
    Gas_Electron_Abundance=Gas_Dict['ElectronAbundance'].astype(np.float64)

    m_p=1.673E-24
    X_H=0.76
    unit_switching=1E10
    mean_molecular_weight=4*m_p/(1+3*X_H+4*X_H*Gas_Electron_Abundance)
    k_B=1.38E-16
    gas_cell_temperature_in_Kelvin=2/3*Gas_Internal_Energy/k_B*unit_switching*mean_molecular_weight

    condition = (gas_cell_temperature_in_Kelvin < 10000)
    index = np.where(condition)[0]

    return np.sum(Gas_Masses[index]),np.sum(Gas_Masses)

def process_data(index):
    snapshot = all_tree['SnapNum'][index]
    subfindid = all_tree['SubfindID'][index]
    cold_gas_mass, gas_mass = Cold_Gas_Mass_Calculator(snapshot, subfindid)
    cold_gas_mass_array[snapshot] += cold_gas_mass

if __name__ == "__main__":
    num_processes = 32  # Set the number of processes

    cold_gas_mass_array = multiprocessing.Array('d', 100)  # Shared array for storing results

    with multiprocessing.Pool(processes=num_processes) as pool:
        indices = list(range(len(all_tree['SnapNum'])))
        list(tqdm(pool.imap(process_data, indices), total=len(indices)))

    cold_gas_mass_array = np.array(cold_gas_mass_array)
