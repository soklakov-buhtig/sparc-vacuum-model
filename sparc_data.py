import zipfile
import numpy as np
from sparc_parser import get_galaxy_geometry_from_table

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



# АВТОМАТИЧЕСКАЯ СБОРКА БАЗЫ ДАННЫХ ИЗ ZIP-АРХИВА НА ЛЕТУ
ARCHIVE_PATH = "Rotmod_LTG.zip"
SPARC_DATABASE = {}


# Список галактик, которые мы хотим собрать из zip-архива
# Автоматически извлекаем имена всех галактик, у которых в архиве есть файл _rotmod.dat
with zipfile.ZipFile(ARCHIVE_PATH, 'r') as z:
    file_names = z.namelist()
    target_galaxies = [f.replace('_rotmod.dat', '') for f in file_names if f.endswith('_rotmod.dat')]
    target_galaxies = sorted(list(set(target_galaxies))) # Сортируем по алфавиту


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
