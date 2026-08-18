import zipfile

def read_galaxy_from_zip(archive_path, galaxy_name):
    # Формируем точное имя файла внутри архива, например "NGC3198_rotmod"
    file_name = f"{galaxy_name.replace(' ', '')}_rotmod.dat"
    
    # Открываем zip-архив в режиме чтения
    with zipfile.ZipFile(archive_path, 'r') as z:
        # Читаем сырые байты файла и декодируем их в текст
        with z.open(file_name) as f:
            lines = [line.decode('utf-8') for line in f.readlines()]
            
    # Списки для сбора числовых колонок
    r_list, vobs_list, vgas_list, vdisk_list = [], [], [], []
    
    # Парсим строки (пропускаем заголовки, если они начинаются не с цифр)
    for line in lines:
        parts = line.strip().split()
        if not parts or not parts[0].replace('.', '', 1).isdigit():
            continue  # Пропускаем строки с текстом заголовков
            
        # Заполняем списки (обычно в файле 4 базовые колонки: Rad, Vobs, Vgas, Vdisk)
        r_list.append(float(parts[0]))
        vobs_list.append(float(parts[1]))
        vgas_list.append(float(parts[2]))
        vdisk_list.append(float(parts[3]))
        
    return (np.array(r_list), np.array(vobs_list), 
            np.array(vgas_list), np.array(vdisk_list))
