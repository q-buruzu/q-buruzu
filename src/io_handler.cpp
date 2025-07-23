#include "io_handler.h"
#include "main_utils.h"

#include <atomic>
#include <iostream>
#include <string>

extern std::atomic<bool> running;
std::atomic<bool> measureRequest(false);
std::atomic<bool> pauseRequest(false);
std::atomic<bool> unpauseRequest(false);

void inputListener() {
        std::string command;

	while (true) {
        	std::getline(std::cin, command);

		if (command == "measure") {
			measureRequest = true;
		} else if (command == "pause") {
			pauseRequest = true;
		} else if (command == "unpause") {
			unpauseRequest = true;
		}
	}
}

void inputHandler() {
	if (measureRequest) {
		measure();
		measureRequest = false;
		pauseRequest = true;
	}

	if (pauseRequest) {
		running = false;
		pauseRequest = false;
	}

	if (unpauseRequest) {
		running = true;
		unpauseRequest = false;
	}
}
