import zipfile
import numpy as np

def read_galaxy_from_zip(archive_path, galaxy_name):
    # Точное имя файла внутри архива, как мы увидели в проводнике
    file_name = f"{galaxy_name.replace(' ', '')}_rotmod.dat"
    
    with zipfile.ZipFile(archive_path, 'r') as z:
        with z.open(file_name) as f:
            lines = [line.decode('utf-8') for line in f.readlines()]
            
    r_list, vobs_list, vgas_list, vdisk_list = [], [], [], []
    
    for line in lines:
        parts = line.strip().split()
        # Пропускаем пустые строки или строки заголовков (комментариев)
        if not parts or parts[0].startswith('#') or not parts[0].replace('.', '', 1).isdigit():
            continue
            
        r_list.append(float(parts[0]))
        vobs_list.append(float(parts[1]))
        vgas_list.append(float(parts[3]))
        vdisk_list.append(float(parts[4]))
        
    return np.array(r_list), np.array(vobs_list), np.array(vgas_list), np.array(vdisk_list)

def get_galaxy_geometry_from_table(galaxy_name):
    search_name = galaxy_name.replace(' ', '')
    
    with open("sparc_to_correct.py", "r", encoding="utf-8") as f:
        lines = f.readlines()
        
    for line in lines:
        parts = line.strip().split()
        if not parts:
            continue
            
        if parts[0] == search_name:
            # 1. Считываем тип галактики (2-я колонка в таблице)
            g_type = int(parts[1])
            
            # 2. Точный масштабный радиус диска Rdisk (12-я колонка, индекс 11)
            r_disk = float(parts[11])
            
            # 3. Физический расчет полутолщины z_c в зависимости от морфологии
            if g_type >= 9:  # Карликовые неправильные системы (Dwarf Irr)
                z_c = r_disk * 0.25  # Диск у карликов значительно толще
            elif g_type <= 3:  # Ранние спирали с балджем
                z_c = r_disk * 0.08  # Тонкие плотные диски
            else:
                z_c = r_disk * 0.12  # Стандартные спирали
                
            # 4. Начальное приближение черной дыры
            m_bh = 0.001 if "3198" in search_name else 0.0001
            
            return r_disk, z_c, m_bh
            
    return 4.0, 0.35, 0.001



# АВТОМАТИЧЕСКАЯ СБОРКА БАЗЫ ДАННЫХ ИЗ ZIP-АРХИВА НА ЛЕТУ
ARCHIVE_PATH = "Rotmod_LTG.zip"
SPARC_DATABASE = {}


# Список галактик, которые мы хотим собрать из zip-архива
target_galaxies = ["NGC 3198", "NGC 2403", "NGC 7331", "NGC 2903", "IC 2574", "NGC 2841", "NGC 3992", "NGC 5055", "NGC 5033", "NGC 4157", "NGC 3877", "NGC 3521"]

for g_name in target_galaxies:
    try:
        # Автоматически вытаскиваем R_d (r_disk) и z_c из текстовой таблицы!
        r_disk, z_c, m_bh = get_galaxy_geometry_from_table(g_name)
        
        r, vobs, vgas, vdisk = read_galaxy_from_zip(ARCHIVE_PATH, g_name)
        
        vbar = np.sqrt(vgas**2 + 0.5 * vdisk**2)
        
        SPARC_DATABASE[g_name] = {
            "R": r,
            "Vobs": vobs,
            "Vgas": vgas,
            "Vdisk": vdisk,           
            "R_d": r_disk,  # Передаем динамический параметр из таблицы
            "z_c": z_c,     # Передаем динамический параметр
            "M_bh": m_bh
        }
    except Exception as e:
        print(f"Ошибка чтения данных для {g_name}: {e}")
