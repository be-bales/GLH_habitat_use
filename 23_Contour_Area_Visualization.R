## contour area visualization

areas <- read.csv("Data/updated_contour_areas.csv") # calculated in ArcGIS Pro

head(areas)

areas <- areas %>%
  dplyr::select(-atotal_50, -atotal_75, -atotal_90)

areas_long <- areas %>%
  pivot_longer(
    starts_with("p"),
    names_to = "Contour",
    values_to = "Percent") %>%
  mutate(
    Contour = factor(Contour,
      levels = c("ptotal_50","ptotal_75","ptotal_90"),
      labels = c("50%","75%","90%")),
   Model = factor(Model, 
                  levels = c("sdmTMB", "MaxEnt")))


cols <- c(
  "Historical" = "black",
  "1-1.9" = "#00BFC4",
  "2-4.5" = "#619CFF",
  "5-8.5" = "#F8766D"
)

shape_vals <- c(
  "Historical" = 16,
  "1-1.9" = 16,
  "2-4.5" = 15,
  "5-8.5" = 17
)

linetype_vals <- c(
  "Historical" = "solid",
  "1-1.9" = "solid",
  "2-4.5" = "dashed",
  "5-8.5" = "dotted"
)

decade_levels <- c("2000-09", "2010-19", "2020-29",
                   "2030-39", "2040-49")

ssp_levels <- c("Historical", "1-1.9", "2-4.5", "5-8.5")

## make connecting line from 2010-2019 to 2020-2029
connectors <- areas_long %>%
  filter(Decade == "2010-19", SSP== "Historical")

con1 <- connectors %>% mutate(SSP = "1-1.9")
con2 <- connectors %>% mutate(SSP = "2-4.5")
con3 <- connectors %>% mutate(SSP = "5-8.5")

areas_long_connected <- bind_rows(areas_long, con1, con2, con3) %>%
  mutate(SSP = factor(SSP, levels = ssp_levels))

ggplot(data = areas_long_connected, aes(x = Decade, y = Percent, color = SSP, group = SSP)) +
  geom_line(aes(linetype = SSP), linewidth = 0.8) +
  geom_point(data = areas_long, aes(shape = SSP), size = 2) +
  
  facet_wrap(~ Model + Contour, nrow = 2,
             labeller = label_wrap_gen(multi_line = FALSE)) +
  scale_color_manual(values = cols, name = "SSP scenario") +
  scale_shape_manual(values = shape_vals, name = "SSP scenario") +
  scale_linetype_manual(values = linetype_vals, name = "SSP scenario") +
  scale_y_continuous(limits = c(0,50)) +
  labs(
    x = "Decade",
    y = "Percent of total occupied area (%)") +
  theme_bw() +
  guides(shape = "none") +
  theme(
    legend.text = element_text(size = 9),
    panel.spacing = unit(0.1, "cm"),
    legend.box.spacing = unit(0.1, "cm"),
    legend.position = "bottom",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    text = element_text(size = 9),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size= 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size= 8),
    strip.text.x = element_text(size =9),
    strip.background = element_rect(fill = "grey90", color = "black")
  )

ggsave("Figures/prediction_contour_change_50_frontiers.tif",
       dpi = 300,
       width = 18,
       height = 10,
       device = "tiff",
       compression = "lzw",
       units = "cm",
       bg = "white")
