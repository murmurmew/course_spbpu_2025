#испр. вар.7 Зайченкова
import math
import json
import matplotlib.pyplot as plt


class WellInjector:
    def __init__(self, params):

        self.gamma_water = params["gamma_water"]  # Относительная плотность воды
        self.md_vdp = params["md_vdp"]  # Глубина до перфорации, м
        self.d_tub = params["d_tub"]  # Диаметр НКТ, м
        self.angle = params["angle"]  # Угол наклона, градусы
        self.roughness = params["roughness"]  # Шероховатость, м
        self.p_wh = params["p_wh"]  # Буферное давление, атм
        self.t_wh = params["t_wh"]  # Температура на устье, °C
        self.temp_grad = params["temp_grad"]  # Геотермический градиент, °C/100м
        self.g = 9.81
        self.rho_water_ref = 1000
        self.pa_to_atm = 101325

    def calculate_water_density(self):
        return self.gamma_water * self.rho_water_ref

    def calculate_temperature(self, depth):
        return self.t_wh + (depth / 100) * self.temp_grad

    def calculate_friction_factor(self, reynolds_number):
        if reynolds_number <= 2000:
            # Ламинарный режим - формула Хагена-Пуазейля
            return 64 / reynolds_number
        else:
            # Турбулентный режим - формула Блазиуса
            return 0.316 / (reynolds_number ** 0.25)

    def calculate_pressure_gradient(self, pressure, temperature, flow_rate, depth_increment=1):
        # Плотность воды
        rho_water = self.calculate_water_density()

        # Угол наклона в радианах
        theta_rad = math.radians(self.angle)

        # Площадь поперечного сечения трубы
        cross_section_area = math.pi * (self.d_tub / 2) ** 2

        if flow_rate > 0:
            # Скорость потока (переводим м³/сут в м³/с)
            flow_velocity = flow_rate / (24 * 3600 * cross_section_area)

            # Динамическая вязкость воды, Па·с
            water_viscosity = 0.001

            # Число Рейнольдса
            reynolds = rho_water * flow_velocity * self.d_tub / water_viscosity

            # Коэффициент трения
            friction_factor = self.calculate_friction_factor(reynolds)
        else:
            flow_velocity = 0
            friction_factor = 0

        # Гравитационная составляющая градиента давления, Па/м
        gravitational_gradient = rho_water * self.g * math.sin(theta_rad)

        # Фрикционная составляющая градиента давления, Па/м
        if flow_rate > 0:
            frictional_gradient = friction_factor * (1 / self.d_tub) * (rho_water * flow_velocity ** 2 / 2)
        else:
            frictional_gradient = 0

        # Суммарный градиент давления
        total_gradient = gravitational_gradient + frictional_gradient

        return total_gradient

    def calculate_bottomhole_pressure(self, flow_rate, num_steps=1000):
        # Расчет для нулевого дебита (только гидростатическое давление)
        if flow_rate == 0:
            water_density = self.calculate_water_density()
            theta_rad = math.radians(self.angle)

            # Гидростатическое давление в Паскалях
            hydrostatic_pressure_pa = water_density * self.g * self.md_vdp * math.sin(theta_rad)

            # Переводим в атмосферы
            hydrostatic_pressure_atm = hydrostatic_pressure_pa / self.pa_to_atm

            # Забойное давление = устьевое + гидростатическое
            bottomhole_pressure = self.p_wh + hydrostatic_pressure_atm

            return bottomhole_pressure

        # Расчет для ненулевых дебитов
        step_size = self.md_vdp / num_steps
        current_pressure = self.p_wh  # Начальное давление на устье, атм

        # Интегрирование по глубине скважины
        for step in range(num_steps):
            # Текущая глубина (середина участка для повышения точности)
            current_depth = (step + 0.5) * step_size

            # Температура на текущей глубине
            current_temperature = self.calculate_temperature(current_depth)

            # Градиент давления на участке, Па/м
            pressure_gradient = self.calculate_pressure_gradient(
                current_pressure, current_temperature, flow_rate, step_size
            )

            # Приращение давления на участке, атм
            pressure_increment = (pressure_gradient * step_size) / self.pa_to_atm

            # Обновление текущего давления
            current_pressure += pressure_increment

        return current_pressure


def calculate_vlp_curve(well_params, flow_rate_range):

    # Создание объекта скважины
    well = WellInjector(well_params)

    # Расчет забойных давлений для каждого дебита
    bottomhole_pressures = []

    print("Расчет VLP кривой:")
    print("Дебит (м³/сут) | Забойное давление (атм)")
    print("-" * 45)

    for flow_rate in flow_rate_range:
        bottomhole_pressure = well.calculate_bottomhole_pressure(flow_rate)
        bottomhole_pressures.append(bottomhole_pressure)
        print(f"{flow_rate:12} | {bottomhole_pressure:20.4f}")

    return flow_rate_range, bottomhole_pressures


def plot_vlp_curve(flow_rates, bottomhole_pressures, save_path='vlp_curve.png'):

    plt.figure(figsize=(12, 8))

    # Построение графика
    plt.plot(flow_rates, bottomhole_pressures, 'b-', linewidth=2.5, marker='o',
             markersize=4, markerfacecolor='red', markeredgecolor='darkred')

    # Настройка внешнего вида графика
    plt.xlabel('Дебит жидкости, м³/сут', fontsize=12, fontweight='bold')
    plt.ylabel('Забойное давление, атм', fontsize=12, fontweight='bold')
    plt.title('VLP кривая нагнетательной скважины\n(Вертикальная производительность)',
              fontsize=14, fontweight='bold', pad=20)

    # Добавление сетки
    plt.grid(True, alpha=0.3, linestyle='--')

    # Настройка осей
    plt.xlim(min(flow_rates) - 10, max(flow_rates) + 10)
    plt.ylim(min(bottomhole_pressures) * 0.98, max(bottomhole_pressures) * 1.02)

    # Добавление аннотаций
    plt.annotate(f'Минимальное давление: {min(bottomhole_pressures):.1f} атм',
                 xy=(flow_rates[0], min(bottomhole_pressures)),
                 xytext=(10, 30), textcoords='offset points',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    plt.annotate(f'Максимальное давление: {max(bottomhole_pressures):.1f} атм',
                 xy=(flow_rates[-1], max(bottomhole_pressures)),
                 xytext=(-100, 30), textcoords='offset points',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.7),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    # Улучшение layout
    plt.tight_layout()

    # Сохранение графика
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"\nГрафик сохранен как: {save_path}")


def save_results_to_json(flow_rates, bottomhole_pressures, filename='output.json'):

    results = {
        "q_liq": flow_rates,
        "p_wf": bottomhole_pressures
    }

    with open(filename, 'w', encoding='utf-8') as file:
        json.dump(results, file, indent=2, ensure_ascii=False)

    print(f"Результаты сохранены в файл: {filename}")


def main():

    # Параметры скважины
    well_parameters = {
        "gamma_water": 0.9867268896045713,  # Относительная плотность воды
        "md_vdp": 2455.488558489585,  # Глубина до перфорации, м
        "d_tub": 0.0674409981763102,  # Диаметр НКТ, м
        "angle": 88.13841355675822,  # Угол наклона, градусы
        "roughness": 0.00020568635171041642,  # Шероховатость, м
        "p_wh": 123.89930690558487,  # Буферное давление, атм
        "t_wh": 22.75053527076139,  # Температура на устье, °C
        "temp_grad": 2.1742238941017424  # Геотермический градиент, °C/100м
    }

    # Диапазон дебитов для анализа
    flow_rates = list(range(0, 401, 10))  # от 0 до 400 м³/сут с шагом 10

    # Расчет VLP кривой
    flow_rates_calculated, bottomhole_pressures = calculate_vlp_curve(well_parameters, flow_rates)

    # Сохранение результатов
    save_results_to_json(flow_rates_calculated, bottomhole_pressures)

    # Построение графика
    plot_vlp_curve(flow_rates_calculated, bottomhole_pressures)

    # Вывод сводной информации
    print("\n" + "=" * 50)
    print("СВОДНАЯ ИНФОРМАЦИЯ:")
    print("=" * 50)
    print(f"Глубина скважины: {well_parameters['md_vdp']:.1f} м")
    print(f"Устьевое давление: {well_parameters['p_wh']:.1f} атм")
    print(f"Диаметр НКТ: {well_parameters['d_tub'] * 1000:.1f} мм")
    print(f"Плотность воды: {well_parameters['gamma_water'] * 1000:.1f} кг/м³")
    print(f"Забойное давление при Q=0: {bottomhole_pressures[0]:.1f} атм")
    print(f"Забойное давление при Q=400: {bottomhole_pressures[-1]:.1f} атм")
    print(f"Прирост давления: {bottomhole_pressures[-1] - bottomhole_pressures[0]:.1f} атм")


# Запуск программы
if __name__ == "__main__":
    main()
