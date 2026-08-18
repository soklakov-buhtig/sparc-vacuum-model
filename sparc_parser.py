def get_galaxy_geometry_from_table(galaxy_name):
    # Приводим имя к верхнему регистру без пробелов для надежности поиска
    search_name = galaxy_name.replace(' ', '').upper()
    
    with open("sparc_to_correct.py", "r", encoding="utf-8") as f:
        lines = f.readlines()
        
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 13: 
            continue
            
        # 1. Проверяем вторую колонку (parts[1]), где записано имя галактики
        if parts[1].upper() == search_name:
            # 2. Считываем тип галактики (3-я колонка, индекс 2)
            g_type = int(parts[2])
            
            # 3. Точный масштабный радиус диска Rdisk (12-я колонка, индекс 11)
            r_disk = float(parts[11])
            
            # 4. Физический расчет полутолщины z_c в зависимости от морфологии
            if g_type >= 9:    # Карликовые неправильные системы (Dwarf Irr)
                z_c = r_disk * 0.25
            elif g_type <= 3:  # Ранние спирали с балджем
                z_c = r_disk * 0.08
            else:              # Стандартные спирали
                z_c = r_disk * 0.12
                
            # 5. Начальное приближение массы центральной черной дыры
            m_bh = 0.001 if "3198" in search_name else 0.0001
            
            return r_disk, z_c, m_bh
            
    # Возврат базовых параметров, если объект не найден в каталоге Table1.mrt
    return 4.0, 0.35, 0.001
