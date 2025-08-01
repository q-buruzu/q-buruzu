#include "random_generator.h"

#include <random>

std::random_device rd;
std::mt19937 gen(rd());

bool randChooseQubit() {
        std::bernoulli_distribution dist(0.5);

        return dist(gen);
}

int randDeviationChoose() {
        std::normal_distribution<> dist(0, 1);

        double value = dist(gen);
        int roundedValue = std::ceil(std::abs(value));
        if (value < 0) roundedValue *= -1;

        return roundedValue;
}

int randChooseGate() {
        std::uniform_int_distribution<> dist(0, 3);

        return dist(gen);
}
