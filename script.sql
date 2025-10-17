create table celestial_bodies
(
    id             bigint unsigned auto_increment
        primary key,
    name           varchar(255) not null,
    mu             double       null,
    radius         double       null,
    semimajor_axis double       null,
    eccentricity   double       null,
    inclination    double       null,
    raan           double       null,
    arg_periapsis  double       null,
    true_anomaly   double       null,
    soi            double       null,
    parent_body    text         null,
    constraint id
        unique (id),
    constraint name
        unique (name)
);


