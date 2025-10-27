import math
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


class WellInjector: #модель скважины
    def __init__(self, params): #параметры из входных данных
        self.gamma_water = params["gamma_water"]
        self.md_vdp = params["md_vdp"]
        self.d_tub = params["d_tub"]
        self.angle = params["angle"]
        self.roughness = params["roughness"]
        self.p_wh = params["p_wh"]
        self.t_wh = params["t_wh"]
        self.temp_grad = params["temp_grad"]

        #физ. константы
        self.g = 9.81
        self.rho_water_ref = 1000
        self.p_atm = 1.01325

    def rho_water(self):
        return self.gamma_water * self.rho_water_ref

    def temperature_profile(self, depth):
        return self.t_wh + (depth / 100) * self.temp_grad

    def pressure_gradient(self, p, t, q_liq, depth_step=1):
        rho = self.rho_water()
        theta_rad = math.radians(self.angle)
        a_tub = math.pi * (self.d_tub / 2) ** 2

        if q_liq > 0:
            v = q_liq / (86400 * a_tub)
            mu_water = 0.001
            re = rho * v * self.d_tub / mu_water

            if re > 0:
                f = self._moody_friction_factor(re)
            else:
                f = 0.04
        else:
            v = 0
            f = 0

        dp_grav = rho * self.g * math.sin(theta_rad) * depth_step
        dp_fric = f * (depth_step / self.d_tub) * (rho * v ** 2 / 2) if q_liq > 0 else 0

        return dp_grav + dp_fric

    def _moody_friction_factor(self, re):
        rel_roughness = self.roughness / self.d_tub
        f = 0.02

        def colebrook_equation(f_guess):
            return 1 / math.sqrt(f_guess) + 2 * math.log10(rel_roughness / 3.7 + 2.51 / (re * math.sqrt(f_guess)))

        try:
            f_solution = fsolve(colebrook_equation, f)[0]
            return max(f_solution, 0.005)
        except:
            return 0.02

    def calculate_bottomhole_pressure(self, q_liq, n_steps=100):
        if q_liq == 0:
            rho = self.rho_water()
            theta_rad = math.radians(self.angle)
            p_hydrostatic = rho * self.g * self.md_vdp * math.sin(theta_rad) / 101325
            return self.p_wh + p_hydrostatic

        depth_step = self.md_vdp / n_steps
        p_current = self.p_wh * 101325

        for i in range(n_steps):
            depth = (i + 0.5) * depth_step
            t_current = self.temperature_profile(depth)
            dp_dz = self.pressure_gradient(p_current / 101325, t_current, q_liq, depth_step)
            delta_p = dp_dz * depth_step
            p_current += delta_p

        return p_current / 101325


def plot_vlp_curve(q_liq_range, p_wf_results):
    """Построение графика VLP кривой"""
    plt.figure(figsize=(10, 6))
    plt.plot(q_liq_range, p_wf_results, 'b-', linewidth=2, marker='o', markersize=4)
    plt.xlabel('Дебит жидкости, м³/сут')
    plt.ylabel('Забойное давление, атм')
    plt.title('VLP кривая нагнетательной скважины')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('vlp_curve.png', dpi=300, bbox_inches='tight')
    plt.show()


def main():
    # Исходные данные
    input_params = {
        "gamma_water": 0.9867268896045713,
        "md_vdp": 2455.488558489585,
        "d_tub": 0.0674409981763102,
        "angle": 88.13841355675822,
        "roughness": 0.00020568635171041642,
        "p_wh": 123.89930690558487,
        "t_wh": 22.75053527076139,
        "temp_grad": 2.1742238941017424
    }

    # Создаем объект скважины
    well = WellInjector(input_params)

    # Диапазон дебитов
    q_liq_range = list(range(0, 400, 10))

    # Расчет забойных давлений
    p_wf_results = []

    print("Расчет VLP для нагнетательной скважины:")
    print("Дебит (м³/сут) | Забойное давление (атм)")
    print("-" * 40)

    for q_liq in q_liq_range:
        p_wf = well.calculate_bottomhole_pressure(q_liq)
        p_wf_results.append(p_wf)
        print(f"{q_liq:12} | {p_wf:20.4f}")

    # Сохранение результатов
    output_data = {
        "q_liq": q_liq_range,
        "p_wf": p_wf_results
    }

    with open("output.json", "w") as f:
        json.dump(output_data, f, indent=2)

    # Построение графика
    plot_vlp_curve(q_liq_range, p_wf_results)
    plt.grid()
    print("\nРезультаты сохранены в файл output.json")
    print("График сохранен в файл vlp_curve.png")


if __name__ == "__main__":
    main()