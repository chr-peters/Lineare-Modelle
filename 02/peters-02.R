# Name: Christian Peters

# No. 3)
# ======

data <- read.csv2('wohnmobil.csv')

regressionFunction <- function(x) {
  0.0005552383 * x - 6.1106419791
}

# a)
plot(data$Verfügbares_Einkommen_pro_Einwohner, data$Wohnmobile_pro_1000_Einwohner,
     xlab = 'Verfuegbares Jahreseinkommen (Euro)', ylab = 'Wohnmobile pro 1000 Einwohner',
     main = 'Wohnmobile nach Einkommen')
curve(regressionFunction, add=TRUE)

# Bewertung:
# Wie im Plot zu sehen, spiegelt die Gerade den Zusammenhang zwischen Jahreseinkommen und Wohnmobilanzahl recht angemessen wieder.
# Ein anderer Zusammenhang als der durch die Gerade repraesentierte lineare Zusammenhang erscheint mir ueberdies wenig plausibel.
# Ich sehe daher keine belastbare Grundlage dazu, die Theorie von Schorsch anzufechten.

# b)
residuals <- by(data, 1:nrow(data), function(row) {
  row$Wohnmobile_pro_1000_Einwohner - regressionFunction(row$Verfügbares_Einkommen_pro_Einwohner)
})

plot(data$Verfügbares_Einkommen_pro_Einwohner, residuals, xlab = 'Verfuegbares Jahreseinkommen (Euro)', ylab = 'Residuum', main = 'Residuen')

# c)
# Get the total number of caravans for each place
caravanCounts <- by(data, 1:nrow(data), function(row) {
  row$Wohnmobile_pro_1000_Einwohner * row$Einwohner / 1000
})
# determine the maximum
maxIndex <- which(caravanCounts == max(caravanCounts))
city <- as.character(data[maxIndex, 'Landkreis_Stadt'])

print(paste0('The most promising place to open the caravan business is ', city, '.'))
# "The most promising place to open the caravan business is Steinfurt."