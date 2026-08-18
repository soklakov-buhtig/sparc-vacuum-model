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
        vgas_list.append(float(parts[2]))
        vdisk_list.append(float(parts[3]))
        
    return np.array(r_list), np.array(vobs_list), np.array(vgas_list), np.array(vdisk_list)

# Физические константы и параметры геометрических масштабов галактик, 
# которые не содержатся внутри rotmod.dat таблиц
GALAXY_PARAMS = {
    "NGC 3198": {"R_d": 4.0, "z_c": 0.35, "M_bh": 0.001},
    "NGC 2403": {"R_d": 2.1, "z_c": 0.25, "M_bh": 0.0001},
    "NGC 7331": {"R_d": 6.2, "z_c": 0.60, "M_bh": 0.090},
    "NGC 2903": {"R_d": 3.0, "z_c": 0.40, "M_bh": 0.015},
    "IC 2574": {"R_d": 2.5, "z_c": 0.20, "M_bh": 0.00001},
    "NGC 2841": {"R_d": 4.1, "z_c": 0.55, "M_bh": 0.120}
}

# АВТОМАТИЧЕСКАЯ СБОРКА БАЗЫ ДАННЫХ ИЗ ZIP-АРХИВА НА ЛЕТУ
ARCHIVE_PATH = "Rotmod_LTG.zip"
SPARC_DATABASE = {}

for g_name, params in GALAXY_PARAMS.items():
    try:
        r, vobs, vgas, vdisk = read_galaxy_from_zip(ARCHIVE_PATH, g_name)
        
        # Рассчитываем суммарный барионный профиль Vbar по канонической формуле SPARC
        # (корень из квадратов газа и диска при эталонном инфракрасном масс-факторе 0.5)
        vbar = np.sqrt(vgas**2 + 0.5 * vdisk**2)
        
        # Рассчитываем барионную массу M_bar в каждой точке (в единицах 10^9 масс Солнца)
        G_const = 4.30091e-6  # pc * (km/s)^2 / M_sun
        m_bar = (vbar**2) * (r * 1000) / G_const / 1e9
        
        # Записываем готовые массивы в итоговый словарь для основного приложения
        SPARC_DATABASE[g_name] = {
            "R": r,
            "Vobs": vobs,
            "Vbar": vbar,
            "M_bar": m_bar,
            "R_d": params["R_d"],
            "z_c": params["z_c"],
            "M_bh": params["M_bh"]
        }
    except Exception as e:
        # Если файл какой-то галактики не нашелся в архиве, сервер выдаст предупреждение, но не упадет
        print(f"Ошибка чтения данных для {g_name}: {e}")
