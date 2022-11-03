download_hpc <- function()
{
    temp = tempfile()
    link = paste("https://archive.ics.uci.edu/ml/machine-learning-databases/",
                 "00235/household_power_consumption.zip", sep = "")
    download.file(link,temp)
    file = "household_power_consumption.txt"
    data = read.table(unz(temp, file), sep = ";", header = TRUE)
    unlink(temp)
    
    data$Index = c(1:nrow(data))
    
    temp_tod = c(1:nrow(data))
    temp_wd = c(1:nrow(data))
    
    for (i in 1:nrow(data)) {
        temp_tod[i] = (i + 1043) %% 1440 + 1
        
        if (i == 1) {
            temp_wd[i] = 6
        } else {
            if (temp_tod[i] < temp_tod[i - 1]) {
                temp_wd[i] = temp_wd[i - 1] %% 7 + 1
            } else {
                temp_wd[i] = temp_wd[i - 1]
            }
        }
    }
    
    data$TimeOfDay = temp_tod
    data$Weekday = temp_wd
    
    gap = data$Global_active_power
    indices = which(gap == "?")
    data = data[-indices, ]
    data$Global_active_power = as.numeric(data$Global_active_power)
    data$Global_reactive_power = as.numeric(data$Global_reactive_power)
    data$Voltage = as.numeric(data$Voltage)
    data$Global_intensity = as.numeric(data$Global_intensity)
    data$Sub_metering_1 = as.numeric(data$Sub_metering_1)
    data$Sub_metering_2 = as.numeric(data$Sub_metering_2)
    
    path = "Numerical Results/Data/UCI/power_consumption.csv"
    write.table(data, path, col.names = TRUE, row.names = FALSE, sep = ",")
}
